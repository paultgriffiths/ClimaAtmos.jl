#=
Tests for the chemistry parameterized tendency.

1. Dispatch: NoChemistry is a no-op
2. Tendency sign: IdealizedChemistry applies a negative tendency to positive tracers
3. Integration: column simulation with chemistry=true steps cleanly, tracers exist and change
=#

using Test
import ClimaComms
ClimaComms.@import_required_backends
import ClimaAtmos as CA
import ClimaCore: Spaces, Fields
using ClimaCore.CommonSpaces
import SciMLBase

include("../test_helpers.jl")

# ============================================================================
# Dispatch Tests
# ============================================================================

@testset "Chemistry Dispatch" begin
    @testset "NoChemistry returns nothing" begin
        result = CA.update_chemistry_sources!(
            nothing, nothing, nothing, nothing,
            CA.NoChemistry(),
        )
        @test isnothing(result)
    end
end

# ============================================================================
# Unit Tests: tendency sign and structure
#
# These tests call update_chemistry_sources! directly on synthetic fields,
# avoiding any dependence on specific τ values or tracer names.
# ============================================================================

@testset "IdealizedChemistry tendency" begin
    FT = Float64
    ᶜspace = ExtrudedCubedSphereSpace(
        FT;
        z_elem = 5,
        z_min = 0,
        z_max = FT(1e4),
        radius = FT(6.4e6),
        h_elem = 4,
        n_quad_points = 4,
        staggering = CellCenter(),
    )

    # Synthetic state: uniform positive tracer values
    ρ    = ones(ᶜspace)
    ρch4 = FT(1e-9) .* ones(ᶜspace)
    ρoh  = FT(1e-12) .* ones(ᶜspace)

    Y  = (; c = (; ρ, ρch4, ρoh))
    Yₜ = (; c = (; ρ = zero(ρ), ρch4 = zero(ρch4), ρoh = zero(ρoh)))

    CA.update_chemistry_sources!(Yₜ, Y, nothing, nothing, CA.IdealizedChemistry())

    @testset "CH4 tendency is negative" begin
        @test all(x -> x < 0, Yₜ.c.ρch4)
    end

    @testset "OH tendency is negative" begin
        @test all(x -> x < 0, Yₜ.c.ρoh)
    end

    @testset "Density tendency is untouched" begin
        @test all(iszero, Yₜ.c.ρ)
    end

    @testset "Tendency is proportional to tracer value" begin
        # dρχ/dt ∝ -ρχ: ratio should be spatially uniform
        rates_ch4 = Yₜ.c.ρch4 ./ Y.c.ρch4
        rates_oh  = Yₜ.c.ρoh  ./ Y.c.ρoh
        @test maximum(rates_ch4) ≈ minimum(rates_ch4)
        @test maximum(rates_oh)  ≈ minimum(rates_oh)
        # Both tracers share the same lifetime
        @test maximum(rates_ch4) ≈ maximum(rates_oh)
    end
end

# ============================================================================
# Integration Tests: structural checks only
# ============================================================================

@testset "IdealizedChemistry Integration" begin
    config = CA.AtmosConfig(
        Dict(
            "config" => "column",
            "initial_condition" => "IsothermalProfile",
            "dt" => "10secs",
            "t_end" => "100secs",
            "z_max" => 30000.0,
            "z_stretch" => false,
            "chemistry" => true,
            "output_default_diagnostics" => false,
        ),
        job_id = "chemistry_column",
    )
    (; Y, p, simulation) = generate_test_simulation(config)

    @testset "Chemistry tracers present in state" begin
        @test :ρch4 ∈ propertynames(Y.c)
        @test :ρoh  ∈ propertynames(Y.c)
    end

    @testset "Chemistry model is IdealizedChemistry" begin
        @test p.atmos.chemistry.chemistry_model isa CA.IdealizedChemistry
    end

    @testset "Tracers change and stay finite after stepping" begin
        ch4_initial = maximum(Y.c.ρch4)
        oh_initial  = maximum(Y.c.ρoh)

        SciMLBase.step!(simulation.integrator, 100.0, true)

        @test maximum(Y.c.ρch4) < ch4_initial   # decaying
        @test maximum(Y.c.ρoh)  < oh_initial
        @test !any(isnan, Y.c.ρch4)
        @test !any(isnan, Y.c.ρoh)
    end
end
