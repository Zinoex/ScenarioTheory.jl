@testset "Sample Discarding" begin
    @testset "Errors" begin
        @test_throws DomainError SampleDiscarding(0, 1, 0)
        @test_throws DomainError SampleDiscarding(-1, 1, 0)
        @test_throws DomainError SampleDiscarding(10, 0, 0)
        @test_throws DomainError SampleDiscarding(10, -1, 0)
        @test_throws DomainError SampleDiscarding(10, 1, -1)
        @test_throws ArgumentError SampleDiscarding(10, 6, 5)
        @test_throws ArgumentError SampleDiscarding(10, 10, 1)
    end

    # 10M samples should be good enough for testing purposes, but we can go higher if needed.
    samples_gen = Data.Integers(1, 10_000)
    decision_vars_gen = Data.Integers(1, 1_000)
    discarded_gen = Data.Integers(0, 1_000)

    sample_decision_discarded_gen = Data.Pairs(samples_gen, Data.Pairs(decision_vars_gen, discarded_gen))

    # Beta in (0, 1], but the exponents of beta are more interesting, so we generate logβ in (-15, 0]
    #  and then exponentiate.
    beta_gen = map(Data.Floats{Float64}(;minimum=-15.0, maximum=0.0, nans=false, infs=false)) do logβ
        return exp(logβ)
    end

    # Check correct range
    sdd_gen = filter(sample_decision_discarded_gen) do (samples, (decision_vars, discarded))
        return decision_vars + discarded <= samples
    end

    Supposition.@check function violation_ranges(sample_decision_discarded=sdd_gen, β=beta_gen)
        samples, (decision_vars, discarded) = sample_decision_discarded
        dist = SampleDiscarding(samples, decision_vars, discarded)
        ϵ = violation(dist, β)

        event!("ϵ", ϵ)

        ϵ[1] == 0.0 && 0.0 <= ϵ[2] <= 1.0
    end

    # Check that if samples == decision_vars, then ϵ = 1.0
    Supposition.@check function samples_equal_decision_vars(samples=samples_gen, β=beta_gen)
        dist = SampleDiscarding(samples, samples, 0)
        ϵ = violation(dist, β)

        event!("ϵ", ϵ)

        ϵ[2] == 1.0
    end

    # More decision variables should lead to higher violation, all else equal.
    sdd_gen = filter(sample_decision_discarded_gen) do (samples, (decision_vars, discarded))
        return decision_vars + discarded < samples
    end

    Supposition.@check function successor_decision_vars(sample_decision_discarded=sdd_gen, β=beta_gen)
        samples, (decision_vars, discarded) = sample_decision_discarded

        dist1 = SampleDiscarding(samples, decision_vars, discarded)
        ϵ1 = violation(dist1, β)
        
        dist2 = SampleDiscarding(samples, decision_vars + 1, discarded)
        ϵ2 = violation(dist2, β)

        event!("ϵ1", ϵ1)
        event!("ϵ2", ϵ2)
        
        ϵ1[2] <= ϵ2[2]
    end

    # More discarded samples should lead to higher violation, all else equal.
    sdd_gen = filter(sample_decision_discarded_gen) do (samples, (decision_vars, discarded))
        return decision_vars + discarded < samples
    end

    Supposition.@check function successor_discarded(sample_decision_discarded=sdd_gen, β=beta_gen)
        samples, (decision_vars, discarded) = sample_decision_discarded

        dist1 = SampleDiscarding(samples, decision_vars, discarded)
        ϵ1 = violation(dist1, β)
        
        dist2 = SampleDiscarding(samples, decision_vars, discarded + 1)
        ϵ2 = violation(dist2, β)

        event!("ϵ1", ϵ1)
        event!("ϵ2", ϵ2)
        
        ϵ1[2] <= ϵ2[2]
    end

    # More samples should lead to lower violation, all else equal.
    sdd_gen = filter(sample_decision_discarded_gen) do (samples, (decision_vars, discarded))
        return decision_vars + discarded <= samples
    end

    Supposition.@check function successor_samples(sample_decision_discarded=sdd_gen, β=beta_gen)
        samples, (decision_vars, discarded) = sample_decision_discarded

        dist1 = SampleDiscarding(samples, decision_vars, discarded)
        ϵ1 = violation(dist1, β)

        dist2 = SampleDiscarding(samples + 1, decision_vars, discarded)
        ϵ2 = violation(dist2, β)

        event!("ϵ1", ϵ1)
        event!("ϵ2", ϵ2)

        ϵ1[2] >= ϵ2[2]
    end

    # Higher β (less confidence) should lead to smaller violation, all else equal.
    beta_pairs = filter(Data.Pairs(beta_gen, beta_gen)) do (β1, β2)
        return β1 < β2
    end

    Supposition.@check function successor_beta(sample_decision_discarded=sdd_gen, β=beta_pairs)
        samples, (decision_vars, discarded) = sample_decision_discarded
        β1, β2 = β

        dist = SampleDiscarding(samples, decision_vars, discarded)
        ϵ1 = violation(dist, β1)
        ϵ2 = violation(dist, β2)

        event!("ϵ1", ϵ1)
        event!("ϵ2", ϵ2)

        ϵ1[2] >= ϵ2[2]
    end

    # Soundness: the true violation should be at most ϵ with confidence at least β.
    Supposition.@check function violation_soundness(sample_decision_discarded=sdd_gen, β=beta_gen)
        samples, (decision_vars, discarded) = sample_decision_discarded
        dist = SampleDiscarding(samples, decision_vars, discarded)
        ϵ = violation(dist, β)[2]

        event!("ϵ", ϵ)

        N = samples
        d = decision_vars
        k = discarded
        r = k + d - 1

        logβ = log(β)
        event!("logβ", logβ)

        # Compute the roundtrip confidence using the formula from the paper.
        logβ_roundtrip = logabsbinomial(r, k)[1] + binomlogcdf(N, ϵ, r)
        event!("logβ_roundtrip", logβ_roundtrip)

        # Check that given a β, we chose an ϵ such that the true violation is at most that much.
        # This corresponds to a higher confidence 1 - β, or β >= β_roundtrip.
        logβ >= logβ_roundtrip
    end
end
