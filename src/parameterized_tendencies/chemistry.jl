function update_chemistry_sources!(_, _, _, _, ::NoChemistry)
    return nothing
end

function update_chemistry_sources!(Yₜ, Y, _, _, ::IdealizedChemistry)
    τ_chem = eltype(Y.c.ρ)(10_000) # ~3 hour lifetime (must be >> dynamical dt ~600s)
    @. Yₜ.c.ρch4 -= Y.c.ρch4 / τ_chem
    @. Yₜ.c.ρoh  -= Y.c.ρoh / τ_chem
    return nothing
end
