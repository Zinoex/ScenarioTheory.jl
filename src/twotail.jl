struct TwoTailViolation <: AbstractViolation
    k::Int
    N::Int

    function TwoTailViolation(k::Int, N::Int)
        if N < 1
            throw(DomainError(N, "expected N ≥ 1"))
        end

        if k < 0
            throw(DomainError(k, "expected k ≥ 0"))
        end
        
        if k > N
            throw(ArgumentError("expected k ≤ N"))
        end
        
        return new(k, N)
    end
end

numfailure(dist::TwoTailViolation) = dist.k
numscenarios(dist::TwoTailViolation) = dist.N

function psi(dist::TwoTailViolation, β::Real, α::Real)
    # $$ \tilde{\Psi}_{k,\delta}(\alpha) = \frac{\delta}{2N} \sum_{m=k}^{N-1} \frac{{m \choose k}}{{N \choose k}} (1-\alpha)^{-(N-m)} + \frac{\delta}{6N} \sum_{m=N+1}^{4N} \frac{{m \choose k}}{{N \choose k}} (1-\alpha)^{m-N} $$
    k = dist.k
    N = dist.N

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

function violation(dist::TwoTailViolation, β::Real; tol=1e-10)
    α_lower = 0
    α_upper = numfailure(dist) / numscenarios(dist)
    while α_upper - α_lower > tol
        α = (α_lower + α_upper) / 2
        left = β / 3 * betainc(numfailure(dist) + 1, numsuccess(dist), α) + β / 6 * betainc(numfailure(dist) + 1, 4 * numscenarios(dist) + 1 - numfailure(dist), α)
        right = (1 + β / (6 * numscenarios(dist))) * α * numscenarios(dist) * (betainc(numfailure(dist), numsuccess(dist) + 1, α) - betainc(numfailure(dist) + 1, numsuccess(dist), α))
        if left > right
            α_upper = α
        else
            α_lower = α
        end
    end
    ϵ_lower = α_lower
    
    if numsuccess(dist) == 0
        ϵ_upper = 1
    else
        α_lower = numfailure(dist) / numscenarios(dist)
        α_upper = 1
        while α_upper - α_lower > tol
            α = (α_lower + α_upper) / 2
            left = (β / 2 - β / 6) * betainc(numfailure(dist) + 1, numsuccess(dist), α) + β / 6 * betainc(numfailure(dist) + 1, 4 * numscenarios(dist) + 1 - numfailure(dist), α)
            right = (1 + β / (6 * numscenarios(dist))) * α * numscenarios(dist) * (betainc(numfailure(dist), numsuccess(dist) + 1, α) - betainc(numfailure(dist) + 1, numsuccess(dist), α))
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

function confidence(dist::TwoTailViolation, ϵ::Tuple{<:Real, <:Real}; tol=1e-10)
    throw(ArgumentError("confidence is not implemented for TwoTailViolation"))
end