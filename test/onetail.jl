@testset "Compression One-tail" begin
    # 10M samples should be good enough for testing purposes, but we can go higher if needed.
    samples_gen = Data.Integers(1, 10_000_000)
    compression_gen = Data.Integers(1, 10_000_000)

    sample_compression_gen = Data.Pairs(samples_gen, compression_gen)

    # Beta in (0, 1], but the exponents of beta are more interesting, so we generate logβ in (-15, 0]
    #  and then exponentiate.
    beta_gen = map(Data.Floats(;minimum=-15.0, maximum=0.0, nans=false, infs=false)) do logβ
        return exp(logβ)
    end

    # Check correct range
    sc_gen = filter(sample_compression_gen) do (samples, compression)
        return compression <= samples
    end

    Supposition.@check function scenario_opt(sample_compression=sc_gen, β=beta_gen)
        samples, compression = sample_compression
        dist = CompressionOneTail(samples, compression)
        ϵ = violation(dist, β)
        ϵ[1] == 0.0 && 0.0 <= ϵ[2] <= 1.0
    end

    # Check that if samples == compression, then ϵ = 1.0
    Supposition.@check function samples_equal_compression(samples=samples_gen, β=beta_gen)
        dist = CompressionOneTail(samples, samples)
        ϵ = violation(dist, β)
        ϵ[2] == 1.0
    end

    # More decision variables should lead to higher violation, all else equal.
    sc_gen = filter(sample_compression_gen) do (samples, compression)
        return compression < samples
    end

    Supposition.@check function successor_compression(sample_compression=sc_gen, β=beta_gen)
        samples, compression = sample_compression

        dist1 = CompressionOneTail(samples, compression)
        ϵ1 = violation(dist1, β)
        
        dist2 = CompressionOneTail(samples, compression + 1)
        ϵ2 = violation(dist2, β)
        
        ϵ1[2] <= ϵ2[2]
    end

    # More samples should lead to lower violation, all else equal.
    sc_gen = filter(sample_compression_gen) do (samples, compression)
        return compression <= samples
    end

    Supposition.@check function successor_samples(sample_compression=sc_gen, β=beta_gen)
        samples, compression = sample_compression

        dist1 = CompressionOneTail(samples, compression)
        ϵ1 = violation(dist1, β)

        dist2 = CompressionOneTail(samples + 1, compression)
        ϵ2 = violation(dist2, β)

        ϵ1[2] >= ϵ2[2]
    end

    # Higher β (less confidence) should lead to smaller violation, all else equal.
    beta_pairs = filter(Data.Pairs(beta_gen, beta_gen)) do (β1, β2)
        return β1 < β2
    end

    Supposition.@check function successor_beta(sample_compression=sc_gen, β=beta_pairs)
        samples, compression = sample_compression
        β1, β2 = β

        dist = CompressionOneTail(samples, compression)
        ϵ1 = violation(dist, β1)
        ϵ2 = violation(dist, β2)

        ϵ1[2] >= ϵ2[2]
    end

    # Soundness: the true violation should be at most ϵ with confidence at least β.
    Supposition.@check function violation_soundness(sample_compression=sc_gen, β=beta_gen)
        samples, compression = sample_compression
        dist = CompressionOneTail(samples, compression)
        ϵ = violation(dist, β)[2]

        # Compute the oracle confidence using the regularized incomplete beta function.
        k = compression
        l = samples - k
        β_roundtrip = ϵ * samples * (ScenarioTheory.betainc(k, l + 1, ϵ) - ScenarioTheory.betainc(k + 1, l, ϵ)) / (1.0 - ScenarioTheory.betainc(l, k + 1, 1.0 - ϵ))
        event!("β_roundtrip", β_roundtrip)

        # Check that given a β, we chose an ϵ such that the true violation is at most that much.
        # This corresponds to a higher confidence 1 - β, or β >= β_roundtrip.
        β >= β_roundtrip
    end
end