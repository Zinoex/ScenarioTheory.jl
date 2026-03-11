@testset "Compression Two-tail" begin
    @testset "Errors" begin
        @test_throws DomainError CompressionTwoTail(0, 0)
        @test_throws ArgumentError CompressionTwoTail(10, 11)
        @test_throws DomainError CompressionTwoTail(-1, 0)
        @test_throws DomainError CompressionTwoTail(10, -1)
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
        dist = CompressionTwoTail(samples, compression)
        ϵ = violation(dist, β)

        event!("ϵ", ϵ)

        0.0 <= ϵ[1] <= 1.0 && 0.0 <= ϵ[2] <= 1.0 && ϵ[1] <= ϵ[2]
    end

    # Check one-tail vs two-tail
    Supposition.@check function one_tail_vs_two_tail(sample_compression=sc_gen, β=beta_gen)
        samples, compression = sample_compression

        dist_one_tail = CompressionOneTail(samples, compression)
        ϵ_one_tail = violation(dist_one_tail, β)

        dist_two_tail = CompressionTwoTail(samples, compression)
        ϵ_two_tail = violation(dist_two_tail, β)

        event!("ϵ_one_tail", ϵ_one_tail)
        event!("ϵ_two_tail", ϵ_two_tail)

        ϵ_one_tail[2] <= ϵ_two_tail[2]
    end

    # Check that if samples == compression, then ϵ = 1.0
    Supposition.@check function samples_equal_compression(samples=samples_gen, β=beta_gen)
        dist = CompressionTwoTail(samples, samples)
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

        dist1 = CompressionTwoTail(samples, compression)
        ϵ1 = violation(dist1, β)
        
        dist2 = CompressionTwoTail(samples, compression + 1)
        ϵ2 = violation(dist2, β)

        event!("ϵ1", ϵ1)
        event!("ϵ2", ϵ2)
        
        ϵ1[1] <= ϵ2[1] && ϵ1[2] <= ϵ2[2]
    end

    # More samples should lead to lower violation, all else equal.
    sc_gen = filter(sample_compression_gen) do (samples, compression)
        return compression <= samples
    end

    Supposition.@check function successor_samples(sample_compression=sc_gen, β=beta_gen)
        samples, compression = sample_compression

        dist1 = CompressionTwoTail(samples, compression)
        ϵ1 = violation(dist1, β)

        dist2 = CompressionTwoTail(samples + 1, compression)
        ϵ2 = violation(dist2, β)

        event!("ϵ1", ϵ1)
        event!("ϵ2", ϵ2)

        ϵ1[1] >= ϵ2[1] && ϵ1[2] >= ϵ2[2]
    end

    # Higher β (less confidence) should lead to smaller violation, all else equal.
    beta_pairs = filter(Data.Pairs(beta_gen, beta_gen)) do (β1, β2)
        return β1 < β2
    end

    Supposition.@check function successor_beta(sample_compression=sc_gen, β=beta_pairs)
        samples, compression = sample_compression
        β1, β2 = β

        dist1 = CompressionTwoTail(samples, compression)
        ϵ1 = violation(dist1, β1)

        dist2 = CompressionTwoTail(samples, compression)
        ϵ2 = violation(dist2, β2)

        event!("ϵ1", ϵ1)
        event!("ϵ2", ϵ2)

        ϵ1[1] <= ϵ2[1] && ϵ2[2] <= ϵ1[2]
    end

    # Soundness: the true violation should fall between ϵ with confidence at least β.
    Supposition.@check function violation_lower_soundness(sample_compression=sc_gen, β=beta_gen)
        samples, compression = sample_compression
        dist = CompressionTwoTail(samples, compression)
        ϵ = violation(dist, β)

        event!("ϵ", ϵ)

        # This will result in NaN β_roundtrip since divisor is zero, but it will only happen if samples == compression == 1.
        if ϵ[1] == 0.0
            reject!()
        end

        N = samples
        k = compression

        # Compute the roundtrip confidence using the regularized incomplete beta function.
        β_lower_roundtrip = ϵ[1] * N * binompdf(N, ϵ[1], k) / (binomccdf(N, ϵ[1], k) / 3 + binomccdf(4 * N + 1 - k, ϵ[1], k) / 6 - ϵ[1] * N * binompdf(N, ϵ[1], k) / (6 * N))
        event!("β_lower_roundtrip", β_lower_roundtrip)

        # Check that given a β, we chose an ϵ such that the true violation is at most that much.
        # This corresponds to a higher confidence 1 - β, or β >= β_lower_roundtrip.
        β >= β_lower_roundtrip
    end

    # Soundness: the true violation should fall between ϵ with confidence at least β.
    Supposition.@check function violation_upper_soundness(sample_compression=sc_gen, β=beta_gen)
        samples, compression = sample_compression
        dist = CompressionTwoTail(samples, compression)
        ϵ = violation(dist, β)

        event!("ϵ", ϵ)

        N = samples
        k = compression

        # Compute the roundtrip confidence using the regularized incomplete beta function.
        β_upper_roundtrip = ϵ[2] * N * binompdf(N, ϵ[2], k) / (binomccdf(N, ϵ[2], k) / 3 + binomccdf(4 * N + 1 - k, ϵ[2], k) / 6 - ϵ[2] * N * binompdf(N, ϵ[2], k) / (6 * N))
        event!("β_upper_roundtrip", β_upper_roundtrip)

        # Check that given a β, we chose an ϵ such that the true violation is at most that much.
        # This corresponds to a higher confidence 1 - β, or β >= β_upper_roundtrip.
        β >= β_upper_roundtrip
    end

    # The true violation should be at most ϵ with confidence approximately 1 - β.
    sc_gen = filter(sample_compression_gen) do (samples, compression)
        return compression < samples
    end

    Supposition.@check function violation_lower_approx(sample_compression=sc_gen, β=beta_gen)
        samples, compression = sample_compression
        dist = CompressionTwoTail(samples, compression)
        ϵ = violation(dist, β)

        event!("ϵ", ϵ)

        # This will result in NaN β_roundtrip since divisor is zero, but it will only happen if samples == compression == 1.
        if ϵ[1] == 0.0
            reject!()
        end

        N = samples
        k = compression

        # Compute the roundtrip confidence using the regularized incomplete beta function.
        β_lower_roundtrip = ϵ[1] * N * binompdf(N, ϵ[1], k) / (binomccdf(N, ϵ[1], k) / 3 + binomccdf(4 * N + 1 - k, ϵ[1], k) / 6 - ϵ[1] * N * binompdf(N, ϵ[1], k) / (6 * N))
        event!("β_lower_roundtrip", β_lower_roundtrip)

        # Check that the computed β_roundtrip is approximately equal to the input β, which would indicate that the violation is tight.
        isapprox(β, β_lower_roundtrip)
    end

    # The true violation should be at most ϵ with confidence approximately 1 - β.
    sc_gen = filter(sample_compression_gen) do (samples, compression)
        return compression < samples
    end

    Supposition.@check function violation_upper_approx(sample_compression=sc_gen, β=beta_gen)
        samples, compression = sample_compression
        dist = CompressionTwoTail(samples, compression)
        ϵ = violation(dist, β)

        event!("ϵ", ϵ)

        N = samples
        k = compression

        # Compute the roundtrip confidence using the regularized incomplete beta function.
        β_upper_roundtrip = ϵ[2] * N * binompdf(N, ϵ[2], k) / (binomccdf(N, ϵ[2], k) / 3 + binomccdf(4 * N + 1 - k, ϵ[2], k) / 6 - ϵ[2] * N * binompdf(N, ϵ[2], k) / (6 * N))
        event!("β_upper_roundtrip", β_upper_roundtrip)

        # Check that the computed β_roundtrip is approximately equal to the input β, which would indicate that the violation is tight.
        isapprox(β, β_upper_roundtrip)
    end
end