function update_chemistry_sources!(_, _, _, _, ::NoChemistry)
    return nothing
end

function update_chemistry_sources!(Yₜ, Y, _, _, ::IdealizedChemistry)
    τ_chem = eltype(Y.c.ρ)(10_000) # ~3 hour lifetime (must be >> dynamical dt ~600s)
    @. Yₜ.c.ρch4 -= Y.c.ρch4 / τ_chem
    @. Yₜ.c.ρoh  -= Y.c.ρoh / τ_chem
    return nothing
end

# ---------------------------------------------------------------------------
# Diurnal sun factor: max(0, cos(lat) × cos(hour_angle))
#   = 1 at local noon on the equator, 0 at night everywhere.
# Dispatches on coordinate type so it works on both sphere and column grids.
# FT is inferred from the type of t (seconds).
# ---------------------------------------------------------------------------
function _sun_factor(coord::Geometry.LatLongZPoint, t::FT) where {FT}
    hour_angle = FT(2π) * mod(t, FT(86400)) / FT(86400) - FT(π)
    return max(FT(0), cosd(coord.lat) * cos(hour_angle))
end
_sun_factor(::Any, t::FT) where {FT} = FT(0.5)  # column/box: time-mean sun

function update_chemistry_sources!(Yₜ, Y, _, t, model::TroposphericChemistry)
    FT      = eltype(Y.c.ρ)
    k_co    = FT(2.4e-13)              # cm³ molec⁻¹ s⁻¹, rate constant CO+OH
    oh_noon = FT(model.oh_noon)        # molecules cm⁻³
    F_CO    = FT(model.co_emission)    # kg m⁻² s⁻¹
    coords  = Fields.coordinate_field(Y.c.ρco)

    # CO loss: dρco/dt -= k_CO × [OH](lat,t) × ρco  (kg m⁻³ s⁻¹)
    # The effective lifetime τ = 1 / (k_CO × oh_noon × sun) varies with lat and time.
    @. Yₜ.c.ρco -= k_co * oh_noon * _sun_factor(coords, FT(t)) * Y.c.ρco

    # Surface CO emission: flux F_CO (kg m⁻² s⁻¹) injected into bottom layer only.
    # Divide by layer thickness (J from local geometry) to get kg m⁻³ s⁻¹.
    # Fields.level is a slice operation, not a pointwise function, so it must be
    # hoisted outside @. (which would otherwise try to broadcast it element-wise).
    ᶜJ = Fields.local_geometry_field(Y.c.ρco).J
    ρco_sfc = Fields.level(Yₜ.c.ρco, 1)
    J_sfc   = Fields.level(ᶜJ, 1)
    @. ρco_sfc += F_CO / J_sfc

    return nothing
end
