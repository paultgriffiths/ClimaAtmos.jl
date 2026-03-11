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

    # Surface CO emission: flux (kg m⁻² s⁻¹) → volumetric tendency (kg m⁻³ s⁻¹).
    # Use DivergenceF2C + SetValue (same pattern as surface_flux_tendency!) to get
    # the correct F/Δz in the bottom cell without dividing by the full 3-D Jacobian.
    if :prescribed_co_emission_timevaryinginput in propertynames(p.tracers)
        TimeVaryingInputs.evaluate!(
            p.tracers.co_emission_prescribed,
            p.tracers.prescribed_co_emission_timevaryinginput,
            t,
        )
        # Wrap scalar emission in C3 field via parent-level copy (same data layout).
        parent(p.tracers.co_emission_C3_sfc) .= parent(p.tracers.co_emission_prescribed)
        btt_co = boundary_tendency_scalar(Y.c.ρco, p.tracers.co_emission_C3_sfc)
        @. Yₜ.c.ρco -= btt_co
    else
        ᶠgradᵥ_co = Operators.GradientC2F()
        ᶜdivᵥ_co  = Operators.DivergenceF2C(
            top    = Operators.SetValue(C3(FT(0))),
            bottom = Operators.SetValue(C3(F_CO)),
        )
        @. Yₜ.c.ρco -= lazy(ᶜdivᵥ_co(FT(0) * ᶠgradᵥ_co(Y.c.ρco)))
    end

    # SO2 + OH → H2SO4  (effective 2nd order, k_SO2 = 9×10⁻¹³ cm³ molec⁻¹ s⁻¹)
    # Stoichiometric mass ratio M_H2SO4 / M_SO2 = 98.08 / 64.06 ≈ 1.531
    k_so2  = FT(9e-13)
    ratio  = FT(98.08 / 64.06)
    @. Yₜ.c.ρso2   -= k_so2 * oh * Y.c.ρso2
    @. Yₜ.c.ρh2so4 += k_so2 * oh * Y.c.ρso2 * ratio

    # Surface SO2 emission (same pattern as CO)
    if :prescribed_so2_emission_timevaryinginput in propertynames(p.tracers)
        TimeVaryingInputs.evaluate!(
            p.tracers.so2_emission_prescribed,
            p.tracers.prescribed_so2_emission_timevaryinginput,
            t,
        )
        parent(p.tracers.so2_emission_C3_sfc) .= parent(p.tracers.so2_emission_prescribed)
        btt_so2 = boundary_tendency_scalar(Y.c.ρso2, p.tracers.so2_emission_C3_sfc)
        @. Yₜ.c.ρso2 -= btt_so2
    end

    return nothing
end
