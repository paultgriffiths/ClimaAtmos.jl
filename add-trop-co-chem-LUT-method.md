# Adding Tropospheric CO Chemistry with a 3-D OH Look-Up Table

## Motivation

Carbon monoxide (CO) is a key trace gas in atmospheric chemistry. Its primary loss pathway in the troposphere is oxidation by the hydroxyl radical (OH):

```
CO + OH → CO₂ + H
```

The rate of this reaction depends on local OH concentration, which varies strongly with latitude, altitude, and season — being highest in the sunlit tropical lower troposphere and near zero at the poles and in the stratosphere.

The goal of this work was to implement a simple but physically meaningful CO chemistry scheme in ClimaAtmos that:

1. Transports a CO tracer through the model atmosphere.
2. Removes CO at each grid point at a rate proportional to the local OH concentration.
3. Sources CO at the surface via a prescribed emission flux.
4. Prescribes OH from a published 3-D climatology (GEOS-Chem / GEOS-5) rather than using a crude approximation.

A previous prototype used a "sun factor" — a simple diurnal scaling of a globally uniform OH value. This was replaced with a **look-up table (LUT)**: a NetCDF file containing monthly-mean OH as a function of longitude, latitude, and altitude, read into the model and interpolated to the current time and grid at each timestep.

---

## Background: How ClimaAtmos Represents Model Options

ClimaAtmos uses Julia's **multiple dispatch** pattern to switch between model variants cleanly. The key idea is:

- An **abstract type** acts as a category label (e.g. `AbstractChemistryModel`).
- **Concrete subtypes** (e.g. `NoChemistry`, `TroposphericChemistry`) are zero-cost tags that carry configuration data.
- Functions are written once per variant. Julia automatically picks the right method based on the type of the argument at compile time — no `if/else` chains needed.

```julia
# Abstract category
abstract type AbstractChemistryModel end

# Concrete variants
struct NoChemistry <: AbstractChemistryModel end
struct IdealizedChemistry <: AbstractChemistryModel end

Base.@kwdef struct TroposphericChemistry <: AbstractChemistryModel
    oh_lut_path::String  = ""       # path to OH NetCDF file
    co_emission::Float64 = 1e-10    # surface flux, kg m⁻² s⁻¹
end
```

`Base.@kwdef` lets you construct the struct with keyword arguments and optional defaults:

```julia
model = TroposphericChemistry(oh_lut_path = "/data/OH_jan.nc", co_emission = 0.0)
```

---

## Files Changed

### 1. `src/solver/types.jl` — Defining the chemistry model type

`TroposphericChemistry` was extended from carrying an `oh_noon` (a single number) to carrying `oh_lut_path` (a path to a NetCDF file) and `co_emission`. These two fields are all the information the model needs to know at startup.

```julia
Base.@kwdef struct TroposphericChemistry <: AbstractChemistryModel
    oh_lut_path::String  = ""       # empty → uniform 1×10⁶ molec cm⁻³ fallback
    co_emission::Float64 = 1e-10    # kg m⁻² s⁻¹
end
```

`AtmosChemistry` (defined just below) wraps this in a named struct so it can be stored alongside water, radiation, and other model components inside the top-level `AtmosModel`.

---

### 2. `config/default_configs/default_config.yml` — Registering the new option

Every configurable option in ClimaAtmos must be declared in `default_config.yml` with a default value. The old `chemistry_oh_noon` key was replaced with:

```yaml
chemistry_oh_lut_path: ""    # path to preprocessed OH NetCDF; empty = uniform fallback
chemistry_co_emission: 1.0e-10
```

An empty string is a safe default — the model will fall back to a spatially uniform OH of 10⁶ molecules cm⁻³, which is a reasonable global average.

---

### 3. `src/solver/type_getters.jl` — Parsing the config into a model object

`type_getters.jl` reads the YAML config and constructs the Julia model objects. The relevant section was updated to read `chemistry_oh_lut_path` instead of `chemistry_oh_noon`:

```julia
TroposphericChemistry(;
    oh_lut_path = parsed_args["chemistry_oh_lut_path"],
    co_emission = parsed_args["chemistry_co_emission"],
)
```

---

### 4. `src/cache/tracer_cache.jl` — Allocating memory for the prescribed OH field

During model initialisation, ClimaAtmos allocates a **cache** — a collection of scratch arrays and pre-loaded data that tendencies can read without re-computing or re-reading files at every timestep.

A new function `oh_cache` was added, following the same pattern as the existing `ozone_cache` function:

```julia
function oh_cache(Y, model::TroposphericChemistry, start_date)
    oh_prescribed = similar(Y.c.ρ)   # allocate a field the same shape as air density

    if isempty(model.oh_lut_path)
        # Fallback: fill with a uniform tropospheric background value
        fill!(oh_prescribed, eltype(oh_prescribed)(1e6))
        return (; oh_prescribed)
    else
        # Load the NetCDF LUT via ClimaUtilities
        extrapolation_bc = (Intp.Periodic(), Intp.Flat(), Intp.Flat())
        prescribed_oh_timevaryinginput = TimeVaryingInput(
            model.oh_lut_path,
            "OH",
            axes(oh_prescribed);
            reference_date = start_date,
            regridder_type = :InterpolationsRegridder,
            regridder_kwargs = (; extrapolation_bc),
            method = LinearInterpolation(),
        )
        return (; oh_prescribed, prescribed_oh_timevaryinginput)
    end
end

# No-op for any other chemistry model (dispatch to this method instead)
oh_cache(_, ::AbstractChemistryModel, _) = (;)
```

**What `TimeVaryingInput` does:** `ClimaUtilities.TimeVaryingInputs.TimeVaryingInput` wraps a NetCDF file and knows how to regrid its contents onto the ClimaCore model grid (the cubed-sphere or box grid) and interpolate in time. `InterpolationsRegridder` performs bilinear interpolation in longitude, latitude, and altitude. `LinearInterpolation()` interpolates between monthly snapshots.

**Extrapolation boundary conditions:** `(Intp.Periodic(), Intp.Flat(), Intp.Flat())` means:
- Longitude: periodic (wraps around at ±180°)
- Latitude: flat (clamps to the nearest edge value beyond ±90°)
- Altitude: flat (clamps rather than extrapolating above or below the LUT range)

The returned named tuple `(; oh_prescribed, prescribed_oh_timevaryinginput)` is merged into the main `p.tracers` cache, making both fields accessible from anywhere in the timestepping loop.

---

### 5. `src/cache/cache.jl` — Wiring the OH cache into model initialisation

A one-line change passes the chemistry model into `tracer_cache` so it can call `oh_cache`:

```julia
tracers = tracer_cache(
    Y, aerosol_names, time_varying_trace_gas_names,
    start_date,
    atmos.chemistry.chemistry_model   # ← new argument
)
```

---

### 6. `src/parameterized_tendencies/chemistry.jl` — The physics: CO loss and emission

This is where the actual atmospheric chemistry happens. The tendency function is called once per timestep and adds or subtracts from the time derivatives of the tracer fields (`Yₜ`).

```julia
function update_chemistry_sources!(Yₜ, Y, p, t, model::TroposphericChemistry)
    FT   = eltype(Y.c.ρ)           # floating-point type (Float32 or Float64)
    k_co = FT(2.4e-13)             # rate constant CO+OH, cm³ molec⁻¹ s⁻¹
    F_CO = FT(model.co_emission)   # surface CO flux, kg m⁻² s⁻¹

    # Step 1: update the OH field from the LUT (if one was loaded)
    if :prescribed_oh_timevaryinginput in propertynames(p.tracers)
        TimeVaryingInputs.evaluate!(
            p.tracers.oh_prescribed,
            p.tracers.prescribed_oh_timevaryinginput,
            t,
        )
    end

    # Step 2: CO loss  dρco/dt -= k_CO × [OH] × ρco
    oh = p.tracers.oh_prescribed   # molecules cm⁻³, interpolated onto model grid
    @. Yₜ.c.ρco -= k_co * oh * Y.c.ρco

    # Step 3: surface CO emission into the bottom model layer only
    ᶜJ      = Fields.local_geometry_field(Y.c.ρco).J
    ρco_sfc = Fields.level(Yₜ.c.ρco, 1)
    J_sfc   = Fields.level(ᶜJ, 1)
    @. ρco_sfc += F_CO / J_sfc

    return nothing
end
```

**Units note:** `ρco` is in kg m⁻³, `[OH]` is in molecules cm⁻³, and `k_CO` is in cm³ molec⁻¹ s⁻¹. The product `k_CO × [OH]` has units s⁻¹ (a loss rate), so `k_CO × [OH] × ρco` has units kg m⁻³ s⁻¹, which is exactly what the tendency `Yₜ.c.ρco` expects.

**`@.` broadcast notation:** The dot syntax applies the operation element-wise across the entire model grid — every column and every level — in a single vectorised call, without writing a loop.

**Surface emission:** `Fields.level(..., 1)` extracts a horizontal slice at the lowest model level. Dividing the flux (kg m⁻² s⁻¹) by the Jacobian `J` (m, the layer thickness) converts it to a volumetric tendency (kg m⁻³ s⁻¹).

---

## The OH Preprocessing Script: `preprocess_oh_lut.jl`

The raw GEOS-Chem OH file (`OH_3Dglobal.geos5.47L.4x5.nc`) has some features that need converting before `ClimaUtilities.InterpolationsRegridder` can consume it:

| Property | Raw file | After preprocessing |
|---|---|---|
| Vertical coordinate | Sigma levels (σ = p/p_surface, 0.99 → 3×10⁻⁵) | Altitude in metres |
| OH units | kg m⁻³ | molecules cm⁻³ |
| Vertical dimension name | `lev` | `z` |

### Vertical coordinate conversion

The raw file's `lev` variable contains **sigma levels** — dimensionless pressure ratios. These decrease monotonically from ~0.99 at the surface to ~3×10⁻⁵ at the model top (~73 km). Converting to altitude uses the isothermal scale-height approximation:

```
z = −H × ln(σ),   H = 7000 m
```

Because sigma decreases upward, the natural log is negative, and the minus sign makes `z` positive and increasing — exactly what the regridder requires.

```julia
H   = 7000.0  # metres, atmospheric scale height
z_m = -H .* log.(lev)   # (47,) array, monotonically increasing
```

### Unit conversion

OH in kg m⁻³ is converted to molecules cm⁻³ using Avogadro's number (Nₐ) and the molar mass of OH (M_OH = 17.008 g mol⁻¹):

```
[molecules cm⁻³] = [kg m⁻³] × (Nₐ / M_OH) × (1 m³ / 10⁶ cm³)
                 ≈ [kg m⁻³] × 3.54 × 10¹⁹
```

### Why no `_FillValue` on the output?

The raw file marks missing data with a `_FillValue`. NCDatasets substitutes `missing` wherever a value equals `_FillValue`. When the preprocessor writes the output, **no `_FillValue` is set**, because zero is a physically valid OH concentration (e.g. in the dark polar stratosphere). Setting `_FillValue = 0.0` would cause the regridder to encounter `missing` values and fail when it tries to convert them to `Float32`.

---

## Test Simulation Config: `config/model_configs/trop_chem_lut.yml`

```yaml
initial_condition: "IsothermalProfile"  # stationary, isothermal atmosphere
start_date: "19850101"                  # matches GEOS-Chem climatology reference date
t_end: "10days"
dt: "600secs"

chemistry: "tropospheric"
chemistry_oh_lut_path: "/path/to/OH_3Dglobal.geos5.47L.4x5_climaatmos.nc"
chemistry_co_emission: 0.0             # pure decay — no source

disable_surface_flux_tendency: true
advection_test: true                   # suppress dynamics; chemistry only
```

Setting `advection_test: true` suppresses all momentum tendencies, leaving pure chemistry. Starting from a spatially uniform CO field, any spatial structure that develops must come entirely from the spatial pattern of OH.

---

## Outcome and Diagnostics

Running for 10 simulated days produces 41 snapshots of `mmrco` (CO mass mixing ratio) at 6-hour intervals. Five diagnostic figures are produced by `plot_co_decay.jl`:

### Figure 1 — Global-mean CO decay
A log-scale time series of domain-averaged CO. The fitted exponential gives an effective e-folding time of **τ ≈ 72 days**, consistent with the global-mean OH burden of the GEOS-Chem climatology. The clean exponential confirms the chemistry scheme is running correctly.

### Figure 2 — Hovmöller diagram
Zonal-mean surface CO normalised by its initial value, plotted as latitude vs time. Tropical latitudes (where OH is highest) show the fastest decay; polar regions (where OH is near zero) show almost no change. This latitudinal structure is entirely driven by the spatial pattern in the OH LUT.

### Figure 3 — Surface CO maps
Three panels: CO at t=0 (spatially uniform), CO at t=10 days (depleted in the tropics), and the ratio CO(t=end)/CO(t=0). The ratio panel directly maps the integrated OH exposure over the 10-day run.

### Figure 4 — Vertical profiles
Domain-mean CO fraction as a function of height, comparing equatorial (|φ|<10°) and polar (|φ|>70°) columns at t=10 days. OH decreases sharply above the tropopause (~12 km) so decay rates fall off at high altitude, visible as the profile bending back toward 1.

### Figure 5 — Surface OH from the LUT
The January OH field from the preprocessed NetCDF plotted at the same size and projection as the ratio panel in Figure 3. Placing them side by side makes the correspondence between high OH and strong CO depletion immediately visible.

---

## How to Reproduce

**Step 1: Preprocess the OH file**

```bash
julia --project preprocess_oh_lut.jl OH_3Dglobal.geos5.47L.4x5.nc
```

**Step 2: Update the LUT path in the config**

Edit `config/model_configs/trop_chem_lut.yml` and set `chemistry_oh_lut_path` to the path of the output file from Step 1.

**Step 3: Run the simulation**

```julia
import ClimaAtmos, SciMLBase
config = ClimaAtmos.AtmosConfig(
    "config/model_configs/trop_chem_lut.yml";
    job_id = "trop_chem_lut_test"
)
sim = ClimaAtmos.AtmosSimulation(config)
SciMLBase.solve!(sim.integrator)
```

Note: `solve!` is not exported from `ClimaAtmos` directly — it must be called via `SciMLBase.solve!`.

**Step 4: Plot**

```bash
julia --project plot_co_decay.jl output/trop_chem_lut_test/output_0001
```

Output figures are written to the current working directory.
