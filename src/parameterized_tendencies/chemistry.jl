import ClimaUtilities.TimeVaryingInputs

function update_chemistry_sources!(_, _, _, _, ::NoChemistry)
    return nothing
end

function update_chemistry_sources!(Yₜ, Y, _, _, ::IdealizedChemistry)
    τ_chem = eltype(Y.c.ρ)(10_000) # ~3 hour lifetime (must be >> dynamical dt ~600s)
    @. Yₜ.c.ρch4 -= Y.c.ρch4 / τ_chem
    @. Yₜ.c.ρoh  -= Y.c.ρoh / τ_chem
    return nothing
end

function update_chemistry_sources!(Yₜ, Y, p, t, model::TroposphericChemistry)
    FT   = eltype(Y.c.ρ)
    k_co = FT(2.4e-13)           # cm³ molec⁻¹ s⁻¹, rate constant CO+OH
    F_CO = FT(model.co_emission) # kg m⁻² s⁻¹

    # Update prescribed OH from time-varying LUT if one was loaded.
    if :prescribed_oh_timevaryinginput in propertynames(p.tracers)
        TimeVaryingInputs.evaluate!(
            p.tracers.oh_prescribed,
            p.tracers.prescribed_oh_timevaryinginput,
            t,
        )
    end

    # CO loss: dρco/dt -= k_CO × [OH](lat,lon,z) × ρco  (kg m⁻³ s⁻¹)
    oh = p.tracers.oh_prescribed  # molecules cm⁻³, on model grid
    @. Yₜ.c.ρco -= k_co * oh * Y.c.ρco

    # Surface CO emission: flux F_CO (kg m⁻² s⁻¹) injected into bottom layer only.
    # Fields.level is a slice, not pointwise, so it must be hoisted outside @.
    ᶜJ      = Fields.local_geometry_field(Y.c.ρco).J
    ρco_sfc = Fields.level(Yₜ.c.ρco, 1)
    J_sfc   = Fields.level(ᶜJ, 1)
    @. ρco_sfc += F_CO / J_sfc

    return nothing
end
