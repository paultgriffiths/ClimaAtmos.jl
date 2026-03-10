#!/usr/bin/env julia
#=
plot_co_decay.jl
================
Produces two figures from a TroposphericChemistry run:

  Figure 1 — Hovmöller: zonal-mean surface-layer CO (kg/kg) vs latitude and time.
             Shows fast equatorial decay and near-constant polar CO.

  Figure 2 — Vertical profile: domain-mean CO mixing ratio vs height at t=0 and t=5d.
             Shows how the surface emission source builds up the bottom layer
             while upper levels decay uniformly.

Run the demo first:

   julia --project=. -e '
     import ClimaAtmos as CA
     config = CA.AtmosConfig("config/model_configs/drybaroclinicwave_chem.yml";
                              job_id = "chem_demo")
     simulation = CA.get_simulation(config)
     CA.solve_atmos!(simulation)
   '

Then plot:

   julia --project=. plot_co_decay.jl output/chem_demo/output_0002

NetCDF variable layout (as written by ClimaAtmos diagnostics):
  mmrco(time, lon, lat, z)   with time in seconds, z in metres
=#

using NCDatasets
using CairoMakie
using Statistics: mean

# ── locate output ────────────────────────────────────────────────────────────
outdir = length(ARGS) > 0 ? ARGS[1] : "output/chem_demo/output_0002"
ncfiles = sort(filter(f -> startswith(basename(f), "mmrco"), readdir(outdir; join=true)))
isempty(ncfiles) && error("No mmrco_*.nc files found in $outdir")

println("Found $(length(ncfiles)) output file(s).")

# ── load data ────────────────────────────────────────────────────────────────
# All time steps live in a single file when output_default_diagnostics writes
# one file per variable.  Shape: (time, lon, lat, z).
NCDataset(ncfiles[1]) do ds
    global t_s, lon, lat, z_m, mmrco
    t_s   = Array(ds["time"])    # seconds since start
    lon   = Array(ds["lon"])
    lat   = Array(ds["lat"])
    z_m   = Array(ds["z"])       # metres
    mmrco = Array(ds["mmrco"])   # (time, lon, lat, z)
end

t_days = t_s ./ 86400
println("Time range: $(round(t_days[1],digits=2)) – $(round(t_days[end],digits=2)) days,  $(length(t_days)) snapshots")
println("Grid: $(length(lon)) lon × $(length(lat)) lat × $(length(z_m)) levels")
println("mmrco at t=0:   min=$(minimum(mmrco[1,:,:,:])) max=$(maximum(mmrco[1,:,:,:]))")
println("mmrco at t=end: min=$(minimum(mmrco[end,:,:,:])) max=$(maximum(mmrco[end,:,:,:]))")

# ── derived quantities ────────────────────────────────────────────────────────
# Zonal mean of the lowest (surface) level  →  (lat, time)
surface_zmean = dropdims(mean(mmrco[:, :, :, 1]; dims=2); dims=2)'   # (lat, time)

# Domain-mean vertical profile at each snapshot  →  (z, time)
profile_mean  = dropdims(mean(mmrco; dims=(2,3)); dims=(2,3))'       # (z, time)

# Normalise by initial value (first snapshot) to show fractional decay
co0_surface = surface_zmean[:, 1]                  # (lat,)
frac_decay  = surface_zmean ./ co0_surface         # (lat, time)

# ── Figure 1: Hovmöller ───────────────────────────────────────────────────────
fig1 = Figure(size = (800, 500))
ax1  = Axis(fig1[1,1];
    xlabel = "Latitude [°]",
    ylabel = "Time [days]",
    title  = "Zonal-mean surface CO / CO₀\n(TroposphericChemistry demo: OH_noon = 5×10⁷ molec cm⁻³, τ_equator ≈ 23 h)",
)

hm = heatmap!(ax1, lat, t_days, frac_decay;
    colormap   = :plasma,
    colorrange = (0, 1),
)
Colorbar(fig1[1,2], hm; label = "CO / CO₀")

save("hovmoller_co_decay.png", fig1)
println("Saved hovmoller_co_decay.png")

# ── Figure 2: Equatorial vs polar vertical profiles at t=5d ──────────────────
# The domain mean hides the lat structure. Instead compare:
#   equatorial column (|lat| < 10°) vs polar column (|lat| > 70°) at t=final,
#   normalised by the initial value.
fig2 = Figure(size = (600, 600))
ax2  = Axis(fig2[1,1];
    xlabel = "CO / CO₀",
    ylabel = "Height [km]",
    title  = "CO depletion profile at t = 5 days\n(equatorial vs polar columns, relative to initial)",
    xscale = identity,
)

z_km = z_m ./ 1000

eq_mask   = abs.(lat) .< 10    # tropical band
pol_mask  = abs.(lat) .> 70    # polar cap

# mmrco shape: (time, lon, lat, z)
co_eq_t0  = mean(mmrco[1,   :, eq_mask,  :]; dims=(1,2))[1,1,:]   # (z,)
co_pol_t0 = mean(mmrco[1,   :, pol_mask, :]; dims=(1,2))[1,1,:]
co_eq_t5  = mean(mmrco[end, :, eq_mask,  :]; dims=(1,2))[1,1,:]
co_pol_t5 = mean(mmrco[end, :, pol_mask, :]; dims=(1,2))[1,1,:]

lines!(ax2, co_eq_t5  ./ co_eq_t0,  z_km; color = :firebrick, linewidth = 2, label = "Equatorial (|φ|<10°)")
lines!(ax2, co_pol_t5 ./ co_pol_t0, z_km; color = :steelblue, linewidth = 2, label = "Polar (|φ|>70°)")
vlines!(ax2, [1.0]; color = :gray, linestyle = :dash, label = "CO₀")
axislegend(ax2; position = :lb)

save("vertical_profile_co.png", fig2)
println("Saved vertical_profile_co.png")

println("\nDone. Open hovmoller_co_decay.png and vertical_profile_co.png.")
