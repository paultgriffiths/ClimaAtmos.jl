"""
preprocess_oh_lut.jl

Converts the GEOS-Chem OH climatology file (OH_3Dglobal.geos5.47L.4x5.nc)
into a format that ClimaUtilities.InterpolationsRegridder can consume:

  1. Vertical coordinate: sigma levels (from the file's `lev` variable) →
     altitude (metres) via the isothermal scale-height formula z = -H × ln(σ),
     H = 7000 m.  The `lev` values are monotonically decreasing sigma (0.99 at
     surface → 3×10⁻⁵ at top), so the resulting z array is monotonically
     increasing as required by ClimaUtilities.
  2. Units: kg m⁻³ → molecules cm⁻³  (factor = Nₐ / M_OH / 1e6)
  3. Dimension names: lon, lat, z, time  (time kept as-is; InterpolationsRegridder
     reads the "hours since ..." attribute directly).

Usage:
    julia --project preprocess_oh_lut.jl /path/to/OH_3Dglobal.geos5.47L.4x5.nc
"""

using NCDatasets

function main(input_path::String, output_path::String)
    # Physical constants
    N_A  = 6.02214076e23  # molecules mol⁻¹
    M_OH = 0.017008       # kg mol⁻¹  (O=16, H=1, → 17.008 g/mol)
    # conversion: [kg m⁻³] → [molecules cm⁻³]
    #   = [kg m⁻³] × (N_A / M_OH) [molecules kg⁻¹] × (1 m³ / 1e6 cm³)
    kg_per_m3_to_molec_per_cm3 = N_A / M_OH / 1.0e6

    println("Reading: $input_path")
    ds_in = NCDataset(input_path, "r")

    lon      = Array(ds_in["lon"])       # (72,)  Float32
    lat      = Array(ds_in["lat"])       # (46,)  Float32
    # Read raw time values (bypass CF DateTime decoding)
    time_raw = ds_in["time"].var[:]      # (12,)  Float64, hours since 1985-01-01
    time_units = ds_in["time"].attrib["units"]
    time_cal   = get(ds_in["time"].attrib, "calendar", "standard")

    # lev contains sigma values (p/p_surf), decreasing from ~0.99 (surface) to ~3e-5 (top).
    lev = Float64.(Array(ds_in["lev"]))  # (47,)

    # NCDatasets returns Union{Missing,Float32} for variables with _FillValue.
    # Read as raw array (bypasses CF fill-value substitution) then handle manually.
    fill_val = Float32(ds_in["OH"].attrib["_FillValue"])
    OH_raw   = ds_in["OH"].var[:, :, :, :]  # (lon=72, lat=46, lev=47, time=12)

    close(ds_in)

    # Altitude from sigma: z = -H × ln(σ), H = 7000 m scale height.
    # lev is ordered surface→top (decreasing σ), so z_m is monotonically increasing.
    H   = 7000.0  # m
    z_m = -H .* log.(lev)
    println("Vertical levels (km): ", round.(z_m ./ 1000, digits=1))

    # Convert units: kg m⁻³ → molecules cm⁻³; replace fill values with 0.
    # Use exact equality to detect fill values (not threshold), so valid small
    # OH concentrations are preserved.
    OH_conv = zeros(Float32, size(OH_raw))
    for idx in eachindex(OH_raw)
        v = OH_raw[idx]
        OH_conv[idx] = (ismissing(v) || v === fill_val) ? 0.0f0 :
                       Float32(Float64(v) * kg_per_m3_to_molec_per_cm3)
    end

    # Write output: dimensions (lon, lat, z, time) — standard CF order
    println("Writing: $output_path")
    ds_out = NCDataset(output_path, "c")

    defDim(ds_out, "lon",  length(lon))
    defDim(ds_out, "lat",  length(lat))
    defDim(ds_out, "z",    length(z_m))
    defDim(ds_out, "time", length(time_raw))  # fixed size — static climatology

    v_lon = defVar(ds_out, "lon",  Float32, ("lon",))
    v_lon.attrib["units"]     = "degrees_east"
    v_lon.attrib["long_name"] = "Longitude"
    v_lon.attrib["axis"]      = "X"
    v_lon[:] = lon

    v_lat = defVar(ds_out, "lat",  Float32, ("lat",))
    v_lat.attrib["units"]     = "degrees_north"
    v_lat.attrib["long_name"] = "Latitude"
    v_lat.attrib["axis"]      = "Y"
    v_lat[:] = lat

    v_z = defVar(ds_out, "z", Float32, ("z",))
    v_z.attrib["units"]     = "m"
    v_z.attrib["long_name"] = "Altitude above sea level"
    v_z.attrib["axis"]      = "Z"
    v_z.attrib["positive"]  = "up"
    v_z[:] = Float32.(z_m)

    # Write raw float time values first, then set CF attributes — avoids NCDatasets
    # CF-variable issue where writing floats to an attrib-decorated variable silently fails.
    v_time = defVar(ds_out, "time", Float64, ("time",))
    v_time[:] = time_raw
    v_time.attrib["units"]    = time_units
    v_time.attrib["calendar"] = time_cal
    v_time.attrib["axis"]     = "T"

    # OH: write as (lon, lat, z, time) — ClimaUtilities expects lon first.
    # NCDatasets reads NetCDF arrays in reversed (Fortran/column-major) order,
    # so OH_conv from Array(ds["OH"]) is already (lon, lat, lev, time).
    # No transposition needed — just copy straight across, renaming lev→z.
    v_OH = defVar(ds_out, "OH", Float32, ("lon", "lat", "z", "time"),
        attrib = ["units" => "molecules cm-3",
                  "long_name" => "Prescribed OH concentration"])
    v_OH[:, :, :, :] = Float32.(OH_conv)

    ds_out.attrib["Title"]       = "GEOS-Chem OH climatology, preprocessed for ClimaAtmos"
    ds_out.attrib["Source"]      = basename(input_path)
    ds_out.attrib["OH_units_original"] = "kg m-3"
    ds_out.attrib["conversion"]  = "molecules cm-3 = kg m-3 * $(kg_per_m3_to_molec_per_cm3)"

    close(ds_out)
    println("Done. Output: $output_path")
    println("OH units: molecules cm⁻³")
    println("Vertical: altitude (m), range $(round(minimum(z_m))) – $(round(maximum(z_m))) m")
end

if length(ARGS) < 1
    error("Usage: julia --project preprocess_oh_lut.jl <input.nc> [output.nc]")
end
input  = ARGS[1]
output = get(ARGS, 2, replace(input, ".nc" => "_climaatmos.nc"))
main(input, output)
