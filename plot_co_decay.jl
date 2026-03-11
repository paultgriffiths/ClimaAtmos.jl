#!/usr/bin/env julia
#=
plot_co_decay.jl
================
Produces four figures from a TroposphericChemistry / OH-LUT run:

  Figure 1 — Global-mean CO vs time (log scale).
              Confirms chemistry is working: exponential decay with the
              characteristic timescale set by the OH climatology.

  Figure 2 — Hovmöller: zonal-mean surface-layer CO fraction vs latitude × time.
              Shows faster tropical decay (high OH) vs slow polar decay (low OH).

  Figure 3 — Surface CO map at t=0 and t=end (lon × lat heatmap).
              At t=0 CO is uniform; at t=end the spatial structure of the OH LUT
              is imprinted.  The ratio panel makes the OH influence explicit.

  Figure 4 — Vertical profile: equatorial vs polar CO fraction at t=end.
              Shows how OH decreases above the tropopause, slowing decay aloft.

Usage:
  julia --project plot_co_decay.jl [output_dir]
  Default output_dir: output/trop_chem_lut_test/output_0001
=#

using NCDatasets
using CairoMakie
using GeoMakie
using GeoMakie.GeoJSON
using NaturalEarth
using Statistics: mean

# ── locate output ────────────────────────────────────────────────────────────
outdir  = length(ARGS) > 0 ? ARGS[1] : "output/trop_chem_lut_test/output_0001"
ncfile  = joinpath(outdir, "mmrco_6h_inst.nc")
isfile(ncfile) || error("File not found: $ncfile")

println("Reading: $ncfile")

# ── load data ────────────────────────────────────────────────────────────────
# Shape: (time, lon, lat, z)
NCDataset(ncfile) do ds
    global t_s, lon, lat, z_m, mmrco
    t_s   = Array(ds["time"])
    lon   = Array(ds["lon"])
    lat   = Array(ds["lat"])
    z_m   = Array(ds["z"])
    mmrco = Array(ds["mmrco"])   # (time, lon, lat, z)
end

t_days = t_s ./ 86400
println("Time range: $(round(t_days[1], digits=2)) – $(round(t_days[end], digits=2)) days,  $(length(t_days)) snapshots")
println("Grid: $(length(lon)) lon × $(length(lat)) lat × $(length(z_m)) levels")

co0   = mmrco[1, :, :, :]   # initial field (lon, lat, z)
co_t0 = mmrco[1,  :, :, 1]  # surface at t=0 (lon, lat)
co_tf = mmrco[end,:, :, 1]  # surface at t=end

# ── Figure 1: Global-mean CO vs time (log scale) ─────────────────────────────
# Domain mean over all levels and columns
co_mean = [mean(mmrco[t, :, :, :]) for t in axes(mmrco, 1)]   # (time,)

fig1 = Figure(size = (700, 420))
ax1  = Axis(fig1[1,1];
    xlabel  = "Time [days]",
    ylabel  = "Global-mean CO mass mixing ratio [kg/kg]",
    title   = "CO decay driven by 3-D OH climatology (GEOS-Chem)",
    yscale  = log10,
)
lines!(ax1, t_days, co_mean; color = :firebrick, linewidth = 2)
scatter!(ax1, t_days, co_mean; color = :firebrick, markersize = 5)

# Overlay a reference exponential fit through the first and last point
τ_fit  = -t_days[end] * 86400 / log(co_mean[end] / co_mean[1])  # seconds
t_ref  = LinRange(t_days[1], t_days[end], 200)
co_ref = co_mean[1] .* exp.(-t_ref .* 86400 ./ τ_fit)
lines!(ax1, t_ref, co_ref; color = :gray, linestyle = :dash,
       label = "Exponential fit  τ ≈ $(round(τ_fit/86400, digits=1)) days")
axislegend(ax1; position = :rt)

save("co_global_mean.png", fig1)
println("Saved co_global_mean.png  (τ_fit ≈ $(round(τ_fit/86400, digits=1)) days)")

# ── Figure 2: Hovmöller ───────────────────────────────────────────────────────
# Zonal-mean surface CO normalised by initial value
surface_zmean = dropdims(mean(mmrco[:, :, :, 1]; dims=2); dims=2)  # (time, lat)
co0_zmean     = surface_zmean[1, :]                                 # (lat,)
frac_hovmoller = surface_zmean ./ co0_zmean'                        # (time, lat)

fig2 = Figure(size = (800, 460))
ax2  = Axis(fig2[1,1];
    xlabel = "Latitude [°]",
    ylabel = "Time [days]",
    title  = "Zonal-mean surface CO / CO₀  (faster tropical decay from high OH)",
)
hm2 = heatmap!(ax2, lat, t_days, frac_hovmoller';
    colormap   = :plasma,
    colorrange = (0.0, 1.0),
)
Colorbar(fig2[1,2], hm2; label = "CO / CO₀")
save("hovmoller_co_decay.png", fig2)
println("Saved hovmoller_co_decay.png")

# ── Figure 3: Surface CO map at t=0 and t=end, plus ratio (with coastlines) ──
coastlines = GeoMakie.coastlines()   # GeoJSON FeatureCollection

function geo_panel!(fig, pos, data_lon, data_lat, data_z; title="", colormap=:YlOrRd_9, colorrange=(0,1))
    ax = GeoAxis(fig[pos...];
        title  = title,
        dest   = "+proj=longlat",
        limits = (-180, 180, -90, 90),
    )
    hm = surface!(ax, data_lon, data_lat, data_z;
        shading        = NoShading,
        colormap       = colormap,
        colorrange     = colorrange,
    )
    lines!(ax, coastlines; color = :black, linewidth = 0.6)
    return ax, hm
end

fig3    = Figure(size = (1000, 820))
co_clim = (0.0, Float64(maximum(co_t0)))
ratio_sfc = co_tf ./ max.(co_t0, eps(Float32))

ax3a, hm3a = geo_panel!(fig3, (1,1), lon, lat, co_t0;
    title      = "Surface CO [kg/kg]  t = $(round(t_days[1],   digits=1)) days",
    colormap   = :YlOrRd_9,
    colorrange = co_clim)

ax3b, hm3b = geo_panel!(fig3, (1,3), lon, lat, co_tf;
    title      = "Surface CO [kg/kg]  t = $(round(t_days[end], digits=1)) days",
    colormap   = :YlOrRd_9,
    colorrange = co_clim)

ax3c, hm3c = geo_panel!(fig3, (2,1:3), lon, lat, ratio_sfc;
    title      = "CO(t=end) / CO(t=0)  — imprint of OH spatial structure",
    colormap   = :RdYlBu_11,
    colorrange = (0.0, 1.0))

Colorbar(fig3[1,2], hm3a; label = "CO [kg/kg]")
Colorbar(fig3[2,4], hm3c; label = "CO / CO₀")

save("surface_co_map.png", fig3)
println("Saved surface_co_map.png")

# ── Figure 4: Vertical profiles at t=end ─────────────────────────────────────
eq_mask  = abs.(lat) .< 10
pol_mask = abs.(lat) .> 70
z_km     = z_m ./ 1000

# mmrco: (time, lon, lat, z)
profile(tidx, mask) = dropdims(mean(mmrco[tidx, :, mask, :]; dims=(1,2)); dims=(1,2))

nt = size(mmrco, 1)
co_eq_t0  = profile(1,  eq_mask)
co_pol_t0 = profile(1,  pol_mask)
co_eq_tf  = profile(nt, eq_mask)
co_pol_tf = profile(nt, pol_mask)

fig4 = Figure(size = (600, 620))
ax4  = Axis(fig4[1,1];
    xlabel = "CO / CO₀",
    ylabel = "Height [km]",
    title  = "CO depletion profile at t = $(round(t_days[end], digits=1)) days\n(equatorial vs polar columns)",
)
lines!(ax4, co_eq_tf  ./ co_eq_t0,  z_km; color = :firebrick,  linewidth = 2, label = "Equatorial (|φ|<10°)")
lines!(ax4, co_pol_tf ./ co_pol_t0, z_km; color = :steelblue,  linewidth = 2, label = "Polar (|φ|>70°)")
vlines!(ax4, [1.0]; color = :gray, linestyle = :dash, label = "No change")
axislegend(ax4; position = :lb)

save("vertical_profile_co.png", fig4)
println("Saved vertical_profile_co.png")

# ── Figure 5: Surface OH from LUT (January, z=1) — same layout as ratio panel ─
# Reads directly from the preprocessed OH LUT so the grid matches the source data
# (lon=72, lat=46 at native GEOS-5 resolution).  The simulation starts 1985-01-01
# so time index 1 (January) is the right slice.
oh_lut_path = joinpath(dirname(outdir), "..", "OH_3Dglobal.geos5.47L.4x5_climaatmos.nc")
oh_lut_path = normpath(oh_lut_path)
if !isfile(oh_lut_path)
    # fallback: look next to the output dir
    oh_lut_path = "OH_3Dglobal.geos5.47L.4x5_climaatmos.nc"
end

NCDataset(oh_lut_path) do ds
    global oh_lon, oh_lat, oh_sfc_jan
    oh_lon      = Array(ds["lon"])          # (72,)
    oh_lat      = Array(ds["lat"])          # (46,)
    oh_sfc_jan  = Array(ds["OH"])[:, :, 1, 1]  # (lon=72, lat=46), surface, January
end

# Figure size chosen to match the width×height of the ratio panel in surface_co_map.png.
# That panel occupies ~1000 × (820/2) px, so use (1000, 460).
fig5 = Figure(size = (1000, 500))
ax5  = GeoAxis(fig5[1,1];
    title  = "Prescribed surface OH [molecules cm⁻³]  (GEOS-Chem climatology, January)",
    dest   = "+proj=longlat",
    limits = (-180, 180, -90, 90),
)
hm5 = surface!(ax5, oh_lon, oh_lat, oh_sfc_jan;
    shading    = NoShading,
    colormap   = :viridis,
    colorrange = (0.0, Float64(maximum(oh_sfc_jan))),
)
lines!(ax5, coastlines; color = :white, linewidth = 0.7)
Colorbar(fig5[1,2], hm5; label = "OH [molecules cm⁻³]")

save("surface_oh_lut.png", fig5)
println("Saved surface_oh_lut.png")

println("\nAll figures saved.")
println("Figures: co_global_mean.png  hovmoller_co_decay.png  surface_co_map.png  vertical_profile_co.png  surface_oh_lut.png")
