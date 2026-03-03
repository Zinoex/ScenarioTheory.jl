struct CompressionOneTail <: AbstractScenarioTheory
    samples::Int
    compressed::Int

    function CompressionOneTail(samples::Int, compressed::Int)
        if samples < 1
            throw(DomainError(samples, "expected N ≥ 1"))
        end

        if compressed < 0
            throw(DomainError(compressed, "expected compressed ≥ 0"))
        end
        
        if compressed > samples
            throw(ArgumentError("expected compressed ≤ samples"))
        end
        
        return new(samples, compressed)
    end
end

function violation(dist::CompressionOneTail, β::Real; tol=1e-10)
    N = dist.samples
    k = dist.compressed

    if N == k
        ϵ = 1.0
    else
        α_lower = 0.0
        α_upper = 1.0
        while α_upper - α_lower > tol
            α = (α_lower + α_upper) / 2

            left = β * binomccdf(N, α, k)
            right = α * N * binompdf(N, α, k)
            
            if left > right
                α_upper = α
            else
                α_lower = α
            end
        end
        ϵ = α_upper
    end
    return zero(ϵ), ϵ
end
