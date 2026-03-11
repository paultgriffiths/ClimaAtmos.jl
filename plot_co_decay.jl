#!/usr/bin/env julia
#=
plot_co_decay.jl
================
Produces eight figures from a TroposphericChemistry / OH-LUT run:

  Figure 1 — Global-mean CO vs time (linear scale).
              Confirms emissions are working: monotonic build-up from zero
              (rate limited by OH chemical loss).

  Figure 2 — Hovmöller: zonal-mean surface-layer CO vs latitude × time.
              Shows spatial structure of the build-up; OH loss moderates tropics.

  Figure 3 — Surface CO map at t≈25% and t=end (lon × lat heatmap with coastlines).

  Figure 4 — Surface CO at t=end — same size and layout as Figure 6 (OH map)
              for direct source/sink/result comparison.

  Figure 5 — Vertical profile: equatorial vs polar CO fraction at t=end.
              Shows how OH decreases above the tropopause, slowing decay aloft.

  Figure 6 — Surface OH from the GEOS-Chem LUT (January, lowest model level).
              Matches Figure 4 exactly in size and layout for direct comparison.

  Figure 7 — Surface CO emissions (CEDS anthropogenic total, January).
              Shows where CO is being emitted at the surface.

  Figure 8 — 3-panel source / sink / result:
              CO emissions | surface OH | CO(t=end)/CO(t=0).
              Designed for direct visual comparison of the drivers and outcome.

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
outdir  = length(ARGS) > 0 ? ARGS[1] : "output/trop_chem_lut_test/output_0007"
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

# ── Figure 1: Global-mean CO vs time ─────────────────────────────────────────
# Domain mean over all levels and columns.
# CO starts at zero and spins up from emissions, so use linear y-axis.
co_mean = [mean(mmrco[t, :, :, :]) for t in axes(mmrco, 1)]   # (time,)

fig1 = Figure(size = (700, 420))
ax1  = Axis(fig1[1,1];
    xlabel  = "Time [days]",
    ylabel  = "Global-mean CO mass mixing ratio [kg/kg]",
    title   = "CO spin-up driven by CEDS emissions + OH loss (GEOS-Chem)",
)
lines!(ax1, t_days, co_mean; color = :firebrick, linewidth = 2)
scatter!(ax1, t_days, co_mean; color = :firebrick, markersize = 5)

save("co_global_mean.png", fig1)
println("Saved co_global_mean.png")

# ── Figure 2: Hovmöller ───────────────────────────────────────────────────────
# Zonal-mean surface CO (absolute) vs latitude × time.
# CO starts at zero; slower polar build-up relative to tropics reflects
# lower OH (less removal) vs. lower emissions at high latitudes.
surface_zmean = dropdims(mean(mmrco[:, :, :, 1]; dims=2); dims=2)  # (time, lat)
co_hov_max    = max(Float64(maximum(surface_zmean)), 1e-20)

fig2 = Figure(size = (800, 460))
ax2  = Axis(fig2[1,1];
    xlabel = "Latitude [°]",
    ylabel = "Time [days]",
    title  = "Zonal-mean surface CO [kg/kg]  (spin-up from zero; CEDS emissions + OH loss)",
)
hm2 = heatmap!(ax2, lat, t_days, surface_zmean';
    colormap   = :plasma,
    colorrange = (0.0, co_hov_max),
)
Colorbar(fig2[1,2], hm2; label = "CO [kg/kg]")
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
    hm = heatmap!(ax, data_lon, data_lat, data_z;
        colormap   = colormap,
        colorrange = colorrange,
    )
    lines!(ax, coastlines; color = :black, linewidth = 0.6)
    return ax, hm
end

# CO starts at zero — base colorrange on the final-time field.
co_clim = (0.0, max(Float64(maximum(co_tf)), 1e-20))

# Figure 3: early time and t=end side by side.
# Use the same colorscale (max at t=end) for both panels so the build-up is visible.
co_t_early = mmrco[max(1, size(mmrco,1)÷4), :, :, 1]   # ~25% through run
fig3 = Figure(size = (1000, 420))

ax3a, hm3a = geo_panel!(fig3, (1,1), lon, lat, co_t_early;
    title      = "Surface CO [kg/kg]  t ≈ $(round(t_days[max(1,end÷4)], digits=1)) days",
    colormap   = :YlOrRd_9,
    colorrange = co_clim)

ax3b, hm3b = geo_panel!(fig3, (1,3), lon, lat, co_tf;
    title      = "Surface CO [kg/kg]  t = $(round(t_days[end], digits=1)) days",
    colormap   = :YlOrRd_9,
    colorrange = co_clim)

Colorbar(fig3[1,2], hm3a; label = "CO [kg/kg]")

save("surface_co_map.png", fig3)
println("Saved surface_co_map.png")

# Figure 4: surface CO at t=end — same layout as Figure 5/6 (OH map) for comparison.
# Since CO starts at zero the ratio CO(t=end)/CO(t=0) is undefined; instead
# show the absolute final CO field to compare directly with the OH and emission maps.
fig_ratio = Figure(size = (1000, 500))
ax_ratio  = GeoAxis(fig_ratio[1,1];
    title  = "Surface CO at t=$(round(t_days[end], digits=0)) days [kg/kg]",
    dest   = "+proj=longlat",
    limits = (-180, 180, -90, 90),
)
hm_ratio = heatmap!(ax_ratio, lon, lat, co_tf;
    colormap   = :YlOrRd_9,
    colorrange = co_clim,
)
lines!(ax_ratio, coastlines; color = :black, linewidth = 0.6)
Colorbar(fig_ratio[1,2], hm_ratio; label = "CO [kg/kg]")

save("co_ratio_map.png", fig_ratio)
println("Saved co_ratio_map.png")

# ── Figure 5: Vertical profiles at t=end ─────────────────────────────────────
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
    xlabel = "CO mass mixing ratio [kg/kg]",
    ylabel = "Height [km]",
    title  = "CO vertical profile at t = $(round(t_days[end], digits=1)) days\n(equatorial vs polar columns)",
)
lines!(ax4, co_eq_tf,  z_km; color = :firebrick,  linewidth = 2, label = "Equatorial (|φ|<10°)")
lines!(ax4, co_pol_tf, z_km; color = :steelblue,  linewidth = 2, label = "Polar (|φ|>70°)")
axislegend(ax4; position = :rb)

save("vertical_profile_co.png", fig4)
println("Saved vertical_profile_co.png")

# ── Figure 6: Surface OH from LUT (January, z=1) — same layout as ratio panel ─
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
hm5 = heatmap!(ax5, oh_lon, oh_lat, oh_sfc_jan;
    colormap   = :viridis,
    colorrange = (0.0, Float64(maximum(oh_sfc_jan))),
)
lines!(ax5, coastlines; color = :white, linewidth = 0.7)
Colorbar(fig5[1,2], hm5; label = "OH [molecules cm⁻³]")

save("surface_oh_lut.png", fig5)
println("Saved surface_oh_lut.png")

# ── Figure 7: Surface CO emissions (CEDS, January) ───────────────────────────
em_lut_path = "CO-em-anthro_CMIP_CEDS_2020_climaatmos.nc"
if !isfile(em_lut_path)
    em_lut_path = joinpath(dirname(outdir), "..", "CO-em-anthro_CMIP_CEDS_2020_climaatmos.nc")
    em_lut_path = normpath(em_lut_path)
end

NCDataset(em_lut_path) do ds
    global em_lon, em_lat, em_sfc_jan
    em_lon     = Array(ds["lon"])               # (720,)
    em_lat     = Array(ds["lat"])               # (360,)
    em_sfc_jan = Array(ds["CO_total"])[:, :, 1] # (lon, lat), January
end

# Display in units of 10⁻¹⁰ kg m⁻² s⁻¹ for readability.
# Log scale: span ~4 orders of magnitude; floor at 0.01 (= 10⁻¹² kg m⁻² s⁻¹)
# so ocean zeros clamp to the minimum colour rather than breaking log(0).
em_scale    = 1f10   # multiply by this to get 10⁻¹⁰ kg m⁻² s⁻¹
em_scaled   = em_sfc_jan .* Float32(em_scale)
em_clim     = (0.01, Float64(maximum(em_scaled)))

fig7 = Figure(size = (1000, 500))
ax7  = GeoAxis(fig7[1,1];
    title  = "Anthropogenic CO emissions (CEDS 2020, January)  [×10⁻¹⁰ kg m⁻² s⁻¹]",
    dest   = "+proj=longlat",
    limits = (-180, 180, -90, 90),
)
hm7 = heatmap!(ax7, em_lon, em_lat, em_scaled;
    colormap   = :YlOrBr_9,
    colorscale = log10,
    colorrange = em_clim,
)
lines!(ax7, coastlines; color = :black, linewidth = 0.6)
Colorbar(fig7[1,2], hm7; label = "CO flux [×10⁻¹⁰ kg m⁻² s⁻¹]", scale = log10)

save("surface_co_emissions.png", fig7)
println("Saved surface_co_emissions.png")

# ── Figure 8: 3-panel source / sink / result ─────────────────────────────────
# Columns: CO emissions | surface OH | CO(t=end)/CO(t=0)
# All panels share the same projection and aspect ratio for direct comparison.
fig8 = Figure(size = (1600, 500))

Label(fig8[0, 1:6], "CO chemistry: source  |  sink  |  result after $(round(t_days[end], digits=0)) days";
    fontsize = 16, font = :bold, tellwidth = false)

# Panel a — CO emissions
ax8a = GeoAxis(fig8[1,1];
    title  = "(a) Emissions [×10⁻¹⁰ kg m⁻² s⁻¹]",
    dest   = "+proj=longlat", limits = (-180, 180, -90, 90),
)
hm8a = heatmap!(ax8a, em_lon, em_lat, em_scaled;
    colormap = :YlOrBr_9, colorscale = log10, colorrange = em_clim)
lines!(ax8a, coastlines; color = :black, linewidth = 0.5)
Colorbar(fig8[1,2], hm8a; label = "flux [×10⁻¹⁰ kg m⁻² s⁻¹]", scale = log10, height = Relative(0.85))

# Panel b — surface OH
ax8b = GeoAxis(fig8[1,3];
    title  = "(b) OH [molecules cm⁻³]  (GEOS-Chem, January)",
    dest   = "+proj=longlat", limits = (-180, 180, -90, 90),
)
hm8b = heatmap!(ax8b, oh_lon, oh_lat, oh_sfc_jan;
    colormap = :viridis,
    colorrange = (0.0, Float64(maximum(oh_sfc_jan))))
lines!(ax8b, coastlines; color = :white, linewidth = 0.5)
Colorbar(fig8[1,4], hm8b; label = "OH [molec cm⁻³]", height = Relative(0.85))

# Panel c — surface CO at t=end (net result of emissions + OH removal)
ax8c = GeoAxis(fig8[1,5];
    title  = "(c) Surface CO at t=$(round(t_days[end], digits=0)) days [kg/kg]",
    dest   = "+proj=longlat", limits = (-180, 180, -90, 90),
)
hm8c = heatmap!(ax8c, lon, lat, co_tf;
    colormap = :YlOrRd_9, colorrange = co_clim)
lines!(ax8c, coastlines; color = :black, linewidth = 0.5)
Colorbar(fig8[1,6], hm8c; label = "CO [kg/kg]", height = Relative(0.85))

save("co_chemistry_diagnosis.png", fig8)
println("Saved co_chemistry_diagnosis.png")

# ── SO2 and H2SO4 plots (require a simulation run with SO2 chemistry active) ──
so2_ncfile    = joinpath(outdir, "mmrso2_6h_inst.nc")
h2so4_ncfile  = joinpath(outdir, "mmrh2so4_6h_inst.nc")
so2_em_path   = "SO2-em-anthro_CMIP_CEDS_2020_climaatmos.nc"

if isfile(so2_ncfile) && isfile(h2so4_ncfile)

    # Load SO2 and H2SO4 output
    NCDataset(so2_ncfile) do ds
        global mmrso2 = Array(ds["mmrso2"])   # (time, lon, lat, z)
    end
    NCDataset(h2so4_ncfile) do ds
        global mmrh2so4 = Array(ds["mmrh2so4"])
    end

    so2_t0  = mmrso2[1,   :, :, 1]   # surface at t=0
    so2_tf  = mmrso2[end, :, :, 1]   # surface at t=end
    h2so4_tf = mmrh2so4[end, :, :, 1]

    # ── Figure 9: Multi-species spin-up time series ───────────────────────────
    # Global-mean surface mixing ratio vs time for CO, SO2, H2SO4.
    # Each species on its own axis (different magnitudes).
    co_mean_sfc   = [mean(mmrco[t,   :, :, 1]) for t in axes(mmrco,   1)]
    so2_mean_sfc  = [mean(mmrso2[t,  :, :, 1]) for t in axes(mmrso2,  1)]
    h2so4_mean_sfc = [mean(mmrh2so4[t,:, :, 1]) for t in axes(mmrh2so4,1)]

    # Log scale is justified — spans several orders of magnitude and H2SO4 grows
    # much slower than CO.  Skip t=0 (all species are zero, undefined on log scale)
    # and replace any remaining non-positive values with NaN so they don't plot.
    nz = x -> [v > 0 ? v : NaN for v in x]
    t_plot = t_days[2:end]
    fig9 = Figure(size = (800, 480))
    ax9  = Axis(fig9[1,1];
        xlabel  = "Time [days]",
        ylabel  = "Global-mean surface MMR [kg kg⁻¹]",
        title   = "Spin-up from zero: CO, SO₂, H₂SO₄  (surface layer)",
        yscale  = log10,
    )
    lines!(ax9, t_plot, nz(co_mean_sfc[2:end]);    color = :firebrick,  linewidth = 2, label = "CO")
    lines!(ax9, t_plot, nz(so2_mean_sfc[2:end]);   color = :steelblue,  linewidth = 2, label = "SO₂")
    lines!(ax9, t_plot, nz(h2so4_mean_sfc[2:end]); color = :darkorange,  linewidth = 2, label = "H₂SO₄")
    axislegend(ax9; position = :lt)
    save("spinup_timeseries.png", fig9)
    println("Saved spinup_timeseries.png")

    # ── Figure 10: Surface SO2 at t=0 and t=end ──────────────────────────────
    so2_clim = (0.0, max(Float64(maximum(so2_tf)), 1e-20))

    so2_t_early = mmrso2[max(1, size(mmrso2,1)÷4), :, :, 1]
    fig10 = Figure(size = (1000, 420))
    ax10a, hm10a = geo_panel!(fig10, (1,1), lon, lat, so2_t_early;
        title      = "Surface SO₂ [kg kg⁻¹]  t ≈ $(round(t_days[max(1,end÷4)], digits=1)) days",
        colormap   = :Blues_9,
        colorrange = so2_clim)
    ax10b, hm10b = geo_panel!(fig10, (1,3), lon, lat, so2_tf;
        title      = "Surface SO₂ [kg kg⁻¹]  t = $(round(t_days[end], digits=1)) days",
        colormap   = :Blues_9,
        colorrange = so2_clim)
    Colorbar(fig10[1,2], hm10a; label = "SO₂ [kg kg⁻¹]")
    save("surface_so2_map.png", fig10)
    println("Saved surface_so2_map.png")

    # ── Figure 11: Surface H2SO4 at t=end ────────────────────────────────────
    h2so4_clim = (0.0, max(Float64(maximum(h2so4_tf)), 1e-20))

    fig11 = Figure(size = (1000, 500))
    ax11 = GeoAxis(fig11[1,1];
        title  = "Surface H₂SO₄ [kg kg⁻¹]  t = $(round(t_days[end], digits=1)) days",
        dest   = "+proj=longlat", limits = (-180, 180, -90, 90),
    )
    hm11 = heatmap!(ax11, lon, lat, h2so4_tf;
        colormap = :Purples_9, colorrange = h2so4_clim)
    lines!(ax11, coastlines; color = :black, linewidth = 0.6)
    Colorbar(fig11[1,2], hm11; label = "H₂SO₄ [kg kg⁻¹]")
    save("surface_h2so4_map.png", fig11)
    println("Saved surface_h2so4_map.png")

    # ── Figure 12: SO2 chemistry diagnosis ───────────────────────────────────
    # 3-panel: SO2 emissions | surface OH | SO2(t=end)
    # Mirrors Figure 8 for CO.
    so2_em_scaled = if isfile(so2_em_path)
        local ds_so2em = NCDataset(so2_em_path)
        local raw = Array(ds_so2em["SO2_total"])[:, :, 1] .* Float32(1f10)
        close(ds_so2em)
        raw
    else
        nothing
    end

    fig12 = Figure(size = (1600, 500))
    Label(fig12[0, 1:6], "SO₂ chemistry: source  |  sink  |  result after $(round(t_days[end], digits=0)) days";
        fontsize = 16, font = :bold, tellwidth = false)

    if !isnothing(so2_em_scaled)
        so2_em_clim = (0.01, Float64(maximum(so2_em_scaled)))
        ax12a = GeoAxis(fig12[1,1];
            title = "(a) SO₂ emissions [×10⁻¹⁰ kg m⁻² s⁻¹]",
            dest = "+proj=longlat", limits = (-180, 180, -90, 90))
        hm12a = heatmap!(ax12a, em_lon, em_lat, so2_em_scaled;
            colormap = :YlOrBr_9,
            colorscale = log10, colorrange = so2_em_clim)
        lines!(ax12a, coastlines; color = :black, linewidth = 0.5)
        Colorbar(fig12[1,2], hm12a; label = "flux [×10⁻¹⁰ kg m⁻² s⁻¹]",
            scale = log10, height = Relative(0.85))
    end

    ax12b = GeoAxis(fig12[1,3];
        title = "(b) OH [molecules cm⁻³]  (GEOS-Chem, January)",
        dest = "+proj=longlat", limits = (-180, 180, -90, 90))
    hm12b = heatmap!(ax12b, oh_lon, oh_lat, oh_sfc_jan;
        colormap = :viridis,
        colorrange = (0.0, Float64(maximum(oh_sfc_jan))))
    lines!(ax12b, coastlines; color = :white, linewidth = 0.5)
    Colorbar(fig12[1,4], hm12b; label = "OH [molec cm⁻³]", height = Relative(0.85))

    ax12c = GeoAxis(fig12[1,5];
        title = "(c) SO₂ surface MMR at t=end [kg kg⁻¹]",
        dest = "+proj=longlat", limits = (-180, 180, -90, 90))
    hm12c = heatmap!(ax12c, lon, lat, so2_tf;
        colormap = :Blues_9, colorrange = so2_clim)
    lines!(ax12c, coastlines; color = :black, linewidth = 0.5)
    Colorbar(fig12[1,6], hm12c; label = "SO₂ [kg kg⁻¹]", height = Relative(0.85))

    save("so2_chemistry_diagnosis.png", fig12)
    println("Saved so2_chemistry_diagnosis.png")

    println("\nSO₂/H₂SO₄ figures: spinup_timeseries.png  surface_so2_map.png  surface_h2so4_map.png  so2_chemistry_diagnosis.png")
else
    println("\nSO₂/H₂SO₄ output not found — re-run the simulation with mmrso2/mmrh2so4 diagnostics enabled.")
    println("Expected: $so2_ncfile")
end

println("\nAll figures saved.")
println("CO figures: co_global_mean.png  hovmoller_co_decay.png  surface_co_map.png  co_ratio_map.png  vertical_profile_co.png  surface_oh_lut.png  surface_co_emissions.png  co_chemistry_diagnosis.png")
