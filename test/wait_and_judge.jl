@testset "Wait-and-judge Scenario Optimization" begin
    @testset "Errors" begin
        @test_throws DomainError WaitAndJudgeScenarioOptimization(0, 0)
        @test_throws ArgumentError WaitAndJudgeScenarioOptimization(10, 11)
        @test_throws DomainError WaitAndJudgeScenarioOptimization(-1, 0)
        @test_throws DomainError WaitAndJudgeScenarioOptimization(10, -1)
    end

    # 10M samples should be good enough for testing purposes, but we can go higher if needed.
    samples_gen = Data.Integers(1, 10_000_000)
    support_gen = Data.Integers(1, 10_000_000)

    sample_support_gen = Data.Pairs(samples_gen, support_gen)

    # Beta in (0, 1], but the exponents of beta are more interesting, so we generate logβ in (-15, 0]
    #  and then exponentiate.
    beta_gen = map(Data.Floats{Float64}(;minimum=-15.0, maximum=0.0, nans=false, infs=false)) do logβ
        return exp(logβ)
    end

    # Check correct range
    ss_gen = filter(sample_support_gen) do (samples, support)
        return support <= samples
    end

    Supposition.@check function wait_and_judge_scenario_opt(sample_support=ss_gen, β=beta_gen)
        samples, support = sample_support
        dist = WaitAndJudgeScenarioOptimization(samples, support)
        ϵ = violation(dist, β)
        ϵ[1] == 0.0 && 0.0 <= ϵ[2] <= 1.0
    end

    # Check that if samples == support, then ϵ = 1.0
    Supposition.@check function samples_equal_support(samples=samples_gen, β=beta_gen)
        dist = WaitAndJudgeScenarioOptimization(samples, samples)
        ϵ = violation(dist, β)
        ϵ[2] == 1.0
    end

    # More decision variables should lead to higher violation, all else equal.
    ss_gen = filter(sample_support_gen) do (samples, support)
        return support < samples
    end

    Supposition.@check function successor_support(sample_support=ss_gen, β=beta_gen)
        samples, support = sample_support

        dist1 = WaitAndJudgeScenarioOptimization(samples, support)
        ϵ1 = violation(dist1, β)
        
        dist2 = WaitAndJudgeScenarioOptimization(samples, support + 1)
        ϵ2 = violation(dist2, β)
        
        ϵ1[2] <= ϵ2[2]
    end

    # More samples should lead to lower violation, all else equal.
    ss_gen = filter(sample_support_gen) do (samples, support)
        return support <= samples
    end

    Supposition.@check function successor_samples(sample_support=ss_gen, β=beta_gen)
        samples, support = sample_support

        dist1 = WaitAndJudgeScenarioOptimization(samples, support)
        ϵ1 = violation(dist1, β)

        dist2 = WaitAndJudgeScenarioOptimization(samples + 1, support)
        ϵ2 = violation(dist2, β)

        ϵ1[2] >= ϵ2[2]
    end

    # Higher β (less confidence) should lead to smaller violation, all else equal.
    beta_pairs = filter(Data.Pairs(beta_gen, beta_gen)) do (β1, β2)
        return β1 < β2
    end

    Supposition.@check function successor_beta(sample_support=ss_gen, β=beta_pairs)
        samples, support = sample_support
        β1, β2 = β

        dist = WaitAndJudgeScenarioOptimization(samples, support)
        ϵ1 = violation(dist, β1)
        ϵ2 = violation(dist, β2)

        ϵ1[2] >= ϵ2[2]
    end

    # Soundness: the true violation should be at most ϵ with confidence at least β.
    Supposition.@check function violation_soundness(sample_support=ss_gen, β=beta_gen)
        samples, support = sample_support
        dist = WaitAndJudgeScenarioOptimization(samples, support)
        ϵ = violation(dist, β)[2]

        N = samples
        k = support

        # Compute the roundtrip confidence using the regularized incomplete beta function.
        β_roundtrip = ϵ * (N + 1) * binompdf(N, ϵ, k) / binomccdf(N + 1, ϵ, k)

        # Check that given a β, we chose an ϵ such that the true violation is at most that much.
        # This corresponds to a higher confidence 1 - β, or β >= β_roundtrip.
        β >= β_roundtrip
    end
end