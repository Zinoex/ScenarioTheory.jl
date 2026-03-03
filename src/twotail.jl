struct CompressionTwoTail <: AbstractScenarioProblem
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

function psi(dist::CompressionTwoTail, β::Real, α::Real)
    # $$ \tilde{\Psi}_{k,\delta}(\alpha) = \frac{\delta}{2N} \sum_{m=k}^{N-1} \frac{{m \choose k}}{{N \choose k}} (1-\alpha)^{-(N-m)} + \frac{\delta}{6N} \sum_{m=N+1}^{4N} \frac{{m \choose k}}{{N \choose k}} (1-\alpha)^{m-N} $$
    N = dist.samples
    k = dist.compressed

    if α < 0 || α > 1
        throw(DomainError(α, "expected α ∈ (0, 1)"))
    end

    T = promote_type(typeof(float(β)), typeof(float(α)))
    if β == 0 || k == N
        return zero(T)
    end

    T = promote_type(T, Float64)
    βT = T(β)
    αT = T(α)

    # Closed form:
    # 
    # sum_{m=k}^{N-1} C(m,k) (1-α)^(-(N-m))
    #   = (1-α)^(k-N) * I_α(k+1, N-k) / α^(k+1)
    # sum_{m=N+1}^{4N} C(m,k) (1-α)^(m-N)
    #   = (1-α)^(k-N) * I_α(k+1, 4N+1-k) / α^(k+1)
    ratio1 = T(betainc(k + 1, N - k, αT)) / αT^(k + 1)
    ratio2 = T(betainc(k + 1, 4 * N + 1 - k, αT)) / αT^(k + 1)
    coeff = (one(T) - αT)^(k - N) / T(binomial(N, k))
    return βT / T(N) * coeff * (ratio1 / 2 + ratio2 / 6)
end

function violation(dist::CompressionTwoTail, β::Real; tol=1e-10)
    N = dist.samples
    k = dist.compressed
    l = N - k

    α_lower = 0
    α_upper = k / N
    while α_upper - α_lower > tol
        α = (α_lower + α_upper) / 2
        left = β / 3 * betainc(k + 1, l, α) + β / 6 * betainc(k + 1, 4 * N + 1 - k, α)
        right = (1 + β / (6 * N)) * α * N * (betainc(k, l + 1, α) - betainc(k + 1, l, α))
        if left > right
            α_upper = α
        else
            α_lower = α
        end
    end
    ϵ_lower = α_lower
    
    if l == 0
        ϵ_upper = 1
    else
        α_lower = k / N
        α_upper = 1
        while α_upper - α_lower > tol
            α = (α_lower + α_upper) / 2
            left = (β / 2 - β / 6) * betainc(k + 1, l, α) + β / 6 * betainc(k + 1, 4 * N + 1 - k, α)
            right = (1 + β / (6 * N)) * α * N * (betainc(k, l + 1, α) - betainc(k + 1, l, α))
            if left > right
                α_upper = α
            else
                α_lower = α
            end
        end
        ϵ_upper = α_upper
    end
    
    return ϵ_lower, ϵ_upper
end