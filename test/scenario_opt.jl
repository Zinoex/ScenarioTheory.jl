@testset "Scenario Optimization" begin
    @testset "Errors" begin
        @test_throws DomainError ScenarioOptimization(0, 0)
        @test_throws ArgumentError ScenarioOptimization(10, 11)
        @test_throws DomainError ScenarioOptimization(-1, 0)
        @test_throws DomainError ScenarioOptimization(10, -1)
    end

    # 10M samples should be good enough for testing purposes, but we can go higher if needed.
    samples_gen = Data.Integers(1, 10_000_000)
    decision_vars_gen = Data.Integers(1, 10_000_000)

    sample_decision_vars_gen = Data.Pairs(samples_gen, decision_vars_gen)

    # Beta in (0, 1], but the exponents of beta are more interesting, so we generate logβ in (-15, 0]
    #  and then exponentiate.
    beta_gen = map(Data.Floats{Float64}(;minimum=-15.0, maximum=0.0, nans=false, infs=false)) do logβ
        return exp(logβ)
    end

    # Check correct range
    sdv_gen = filter(sample_decision_vars_gen) do (samples, decision_vars)
        return decision_vars <= samples
    end

    Supposition.@check function violation_ranges(sample_decision_vars=sdv_gen, β=beta_gen)
        samples, decision_vars = sample_decision_vars
        dist = ScenarioOptimization(samples, decision_vars)
        ϵ = violation(dist, β)

        event!("ϵ", ϵ)

        ϵ[1] == 0.0 && 0.0 <= ϵ[2] <= 1.0
    end

    # Check that if samples == decision_vars, then ϵ = 1.0
    Supposition.@check function samples_equal_decision_vars(samples=samples_gen, β=beta_gen)
        dist = ScenarioOptimization(samples, samples)
        ϵ = violation(dist, β)

        event!("ϵ", ϵ)

        ϵ[2] == 1.0
    end

    # More decision variables should lead to higher violation, all else equal.
    sdv_gen = filter(sample_decision_vars_gen) do (samples, decision_vars)
        return decision_vars < samples
    end

    Supposition.@check function successor_decision_vars(sample_decision_vars=sdv_gen, β=beta_gen)
        samples, decision_vars = sample_decision_vars

        dist1 = ScenarioOptimization(samples, decision_vars)
        ϵ1 = violation(dist1, β)
        
        dist2 = ScenarioOptimization(samples, decision_vars + 1)
        ϵ2 = violation(dist2, β)

        event!("ϵ1", ϵ1)
        event!("ϵ2", ϵ2)
        
        ϵ1[2] <= ϵ2[2]
    end

    # More samples should lead to lower violation, all else equal.
    sdv_gen = filter(sample_decision_vars_gen) do (samples, decision_vars)
        return decision_vars <= samples
    end

    Supposition.@check function successor_samples(sample_decision_vars=sdv_gen, β=beta_gen)
        samples, decision_vars = sample_decision_vars

        dist1 = ScenarioOptimization(samples, decision_vars)
        ϵ1 = violation(dist1, β)

        dist2 = ScenarioOptimization(samples + 1, decision_vars)
        ϵ2 = violation(dist2, β)

        event!("ϵ1", ϵ1)
        event!("ϵ2", ϵ2)

        ϵ1[2] >= ϵ2[2]
    end

    # Higher β (less confidence) should lead to smaller violation, all else equal.
    beta_pairs = filter(Data.Pairs(beta_gen, beta_gen)) do (β1, β2)
        return β1 < β2
    end

    Supposition.@check function successor_beta(sample_decision_vars=sdv_gen, β=beta_pairs)
        samples, decision_vars = sample_decision_vars
        β1, β2 = β

        dist = ScenarioOptimization(samples, decision_vars)
        ϵ1 = violation(dist, β1)
        ϵ2 = violation(dist, β2)

        event!("ϵ1", ϵ1)
        event!("ϵ2", ϵ2)

        ϵ1[2] >= ϵ2[2]
    end

    # Soundness: the true violation should be at most ϵ with confidence at least β.
    Supposition.@check function violation_soundness(sample_decision_vars=sdv_gen, β=beta_gen)
        samples, decision_vars = sample_decision_vars
        dist = ScenarioOptimization(samples, decision_vars)
        ϵ = violation(dist, β)[2]

        event!("ϵ", ϵ)

        N = samples
        k = decision_vars - 1

        # Compute the roundtrip confidence using the regularized incomplete beta function.
        β_roundtrip = binomcdf(N, ϵ, k)
        event!("β_roundtrip", β_roundtrip)

        # Check that given a β, we chose an ϵ such that the true violation is at most that much.
        # This corresponds to a higher confidence 1 - β, or β >= β_roundtrip.
        β >= β_roundtrip
    end
end