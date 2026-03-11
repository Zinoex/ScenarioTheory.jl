struct CompressionTwoTail <: AbstractScenarioTheory
    samples::Int
    compressed::Int

    function CompressionTwoTail(samples::Int, compressed::Int)
        if samples < 1
            throw(DomainError(samples, "expected samples ≥ 1"))
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

function violation(theory::CompressionTwoTail, β::Real; tol_steps=20)
    N = theory.samples
    k = theory.compressed

    α_lower = 0.0
    α_upper = k / N
    while α_upper > nextfloat(α_lower)
        α = (α_lower + α_upper) / 2
        left = (β / 3) * binomccdf(N, α, k) + (β / 6) * binomccdf(4 * N + 1 - k, α, k)
        right = (1 + β / (6 * N)) * α * N * binompdf(N, α, k)
        if left > nextfloat(right, tol_steps)
            α_lower = α
        else
            α_upper = α
        end
    end
    ϵ_lower = α_lower
    
    if N == k
        ϵ_upper = 1.0
    else
        α_lower = k / N
        α_upper = 1.0
        while α_upper > nextfloat(α_lower)
            α = (α_lower + α_upper) / 2
            left = (β / 3) * binomccdf(N, α, k) + (β / 6) * binomccdf(4 * N + 1 - k, α, k)
            right = (1 + β / (6 * N)) * α * N * binompdf(N, α, k)
            if left > nextfloat(right, tol_steps)
                α_upper = α
            else
                α_lower = α
            end
        end
        ϵ_upper = α_upper
    end
    
    return ϵ_lower, ϵ_upper
end