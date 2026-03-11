@testset "Compression One-tail" begin
    @testset "Errors" begin
        @test_throws DomainError CompressionOneTail(0, 0)
        @test_throws ArgumentError CompressionOneTail(10, 11)
        @test_throws DomainError CompressionOneTail(-1, 0)
        @test_throws DomainError CompressionOneTail(10, -1)
    end

    # 10M samples should be good enough for testing purposes, but we can go higher if needed.
    samples_gen = Data.Integers(1, 10_000_000)
    compression_gen = Data.Integers(1, 10_000_000)

    sample_compression_gen = Data.Pairs(samples_gen, compression_gen)

    # Beta in (0, 1], but the exponents of beta are more interesting, so we generate logβ in (-15, 0]
    #  and then exponentiate.
    beta_gen = map(Data.Floats{Float64}(;minimum=-15.0, maximum=0.0, nans=false, infs=false)) do logβ
        return exp(logβ)
    end

    # Check correct range
    sc_gen = filter(sample_compression_gen) do (samples, compression)
        return compression <= samples
    end

    Supposition.@check function violation_ranges(sample_compression=sc_gen, β=beta_gen)
        samples, compression = sample_compression
        dist = CompressionOneTail(samples, compression)
        ϵ = violation(dist, β)

        event!("ϵ", ϵ)

        ϵ[1] == 0.0 && 0.0 <= ϵ[2] <= 1.0
    end

    # Check that if samples == compression, then ϵ = 1.0
    Supposition.@check function samples_equal_compression(samples=samples_gen, β=beta_gen)
        dist = CompressionOneTail(samples, samples)
        ϵ = violation(dist, β)

        event!("ϵ", ϵ)

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

        event!("ϵ1", ϵ1)
        event!("ϵ2", ϵ2)
        
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

        event!("ϵ1", ϵ1)
        event!("ϵ2", ϵ2)

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

        event!("ϵ1", ϵ1)
        event!("ϵ2", ϵ2)

        ϵ1[2] >= ϵ2[2]
    end

    # Soundness: the true violation should be at most ϵ with confidence at least β.
    Supposition.@check function violation_soundness(sample_compression=sc_gen, β=beta_gen)
        samples, compression = sample_compression
        dist = CompressionOneTail(samples, compression)
        ϵ = violation(dist, β)[2]

        event!("ϵ", ϵ)

        N = samples
        k = compression

        # Compute the roundtrip confidence using the regularized incomplete beta function.
        β_roundtrip = ϵ * N * binompdf(N, ϵ, k) / binomccdf(N, ϵ, k)
        event!("β_roundtrip", β_roundtrip)

        # Check that given a β, we chose an ϵ such that the true violation is at most that much.
        # This corresponds to a higher confidence 1 - β, or β >= β_roundtrip.
        β >= β_roundtrip
    end

    # The true violation should be at most ϵ with confidence approximately 1 - β.
    sc_gen = filter(sample_compression_gen) do (samples, compression)
        return compression < samples
    end

    Supposition.@check function violation_approx(sample_compression=sc_gen, β=beta_gen)
        samples, compression = sample_compression
        dist = CompressionOneTail(samples, compression)
        ϵ = violation(dist, β)[2]

        event!("ϵ", ϵ)

        N = samples
        k = compression

        # Compute the roundtrip confidence using the regularized incomplete beta function.
        β_roundtrip = ϵ * N * binompdf(N, ϵ, k) / binomccdf(N, ϵ, k)
        event!("β_roundtrip", β_roundtrip)

        # Check that the computed β_roundtrip is approximately equal to the input β, which would indicate that the violation is tight.
        isapprox(β, β_roundtrip)
    end
end