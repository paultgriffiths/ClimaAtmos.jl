"""
preprocess_so2_emissions.jl

Convert the CMIP/CEDS anthropogenic SO2 emissions file to a format compatible
with ClimaUtilities.TimeVaryingInputs / InterpolationsRegridder.

Input:  SO2-em-anthro_CMIP_CEDS_2020.nc
          - dims:    (lon=720, lat=360, time=12), 0.5° resolution
          - sectors: SO2_agr, SO2_ene, SO2_ind, SO2_tra, SO2_rco, SO2_slv, SO2_wst, SO2_shp
          - units:   kg m⁻² s⁻¹
          - time:    "days since 2020-01-01" (monthly, Jan–Dec)
          - _FillValue = 1e20

Output: SO2-em-anthro_CMIP_CEDS_2020_climaatmos.nc
          - dims:    (lon=720, lat=360, time=12)  — no _FillValue on SO2_total
          - variable: SO2_total  [kg m⁻² s⁻¹]   (sum of all sectors; missing→0)
          - time:    relabelled to "hours since 1985-01-01" to match simulation
                     reference date; first-of-month timestamps (Jan 1 = t=0)

Usage:
    julia --project preprocess_so2_emissions.jl [input_file] [output_file]
"""

using NCDatasets
using Dates

# ── paths ─────────────────────────────────────────────────────────────────────
input_file  = length(ARGS) >= 1 ? ARGS[1] : "SO2-em-anthro_CMIP_CEDS_2020.nc"
output_file = length(ARGS) >= 2 ? ARGS[2] : "SO2-em-anthro_CMIP_CEDS_2020_climaatmos.nc"

println("Input:  $input_file")
println("Output: $output_file")

# ── sector variable names in the CEDS file ────────────────────────────────────
SECTORS = ["SO2_agr", "SO2_ene", "SO2_ind", "SO2_tra", "SO2_rco", "SO2_slv", "SO2_wst", "SO2_shp"]

# ── reference date for output time axis (must match simulation start_date) ────
REF_DATE = DateTime(1985, 1, 1)

NCDataset(input_file, "r") do ds_in

    lon  = Array(ds_in["lon"])   # (720,)
    lat  = Array(ds_in["lat"])   # (360,)
    nlon = length(lon)
    nlat = length(lat)
    nt   = 12

    # ── sum all sectors, replacing missing (1e20 fill) with 0.0 ──────────────
    println("Summing $(length(SECTORS)) sectors …")
    so2_total = zeros(Float32, nlon, nlat, nt)   # (lon, lat, time)

    for sector in SECTORS
        raw = Array(ds_in[sector])   # NCDatasets gives (lon, lat, time) with
                                     # missing where value == _FillValue
        for k in 1:nt, j in 1:nlat, i in 1:nlon
            v = raw[i, j, k]
            so2_total[i, j, k] += ismissing(v) ? 0f0 : Float32(v)
        end
        println("  added $sector  (max = $(maximum(skipmissing(raw))))")
    end
    println("SO2_total range: $(minimum(so2_total)) – $(maximum(so2_total)) kg m⁻² s⁻¹")

    # ── build output time axis ────────────────────────────────────────────────
    # Use the 1st of each calendar month so that the first timestamp (Jan 1)
    # coincides exactly with the simulation start date — ClimaUtilities refuses
    # to evaluate before the first data point.
    month_starts = [DateTime(year(REF_DATE), m, 1) for m in 1:12]
    time_hours = Float64[
        (dt - REF_DATE).value / (1000 * 3600) for dt in month_starts
    ]

    # ── write output ──────────────────────────────────────────────────────────
    println("Writing $output_file …")
    NCDataset(output_file, "c") do ds_out

        # dimensions
        defDim(ds_out, "lon",  nlon)
        defDim(ds_out, "lat",  nlat)
        defDim(ds_out, "time", nt)   # fixed size — avoids unlimited-dim write issues

        # coordinate variables
        v_lon = defVar(ds_out, "lon", Float32, ("lon",))
        v_lon[:] = Float32.(lon)
        v_lon.attrib["units"]         = "degrees_east"
        v_lon.attrib["standard_name"] = "longitude"

        v_lat = defVar(ds_out, "lat", Float32, ("lat",))
        v_lat[:] = Float32.(lat)
        v_lat.attrib["units"]         = "degrees_north"
        v_lat.attrib["standard_name"] = "latitude"

        v_time = defVar(ds_out, "time", Float64, ("time",))
        v_time[:] = time_hours
        v_time.attrib["units"]    = "hours since $(Dates.format(REF_DATE, "yyyy-mm-dd HH:MM:SS"))"
        v_time.attrib["calendar"] = "standard"

        # SO2 total — no _FillValue: zero is a valid emission (ocean grid cells)
        v_SO2 = defVar(ds_out, "SO2_total", Float32, ("lon", "lat", "time"))
        v_SO2[:, :, :] = so2_total
        v_SO2.attrib["units"]     = "kg m-2 s-1"
        v_SO2.attrib["long_name"] = "Total anthropogenic SO2 emission (all CEDS sectors)"
    end

    println("Done. Global total (Jan): $(sum(so2_total[:,:,1]) * 0.5^2 * (π/180)^2 * 6.371e6^2) kg s⁻¹")
end
