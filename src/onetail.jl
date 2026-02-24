struct OneTailViolation <: AbstractViolation
    k::Int
    N::Int

    function OneTailViolation(k::Int, N::Int)
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

numfailure(dist::OneTailViolation) = dist.k
numscenarios(dist::OneTailViolation) = dist.N

function psi(dist::OneTailViolation, β::Real, α::Real)
    # $$ \Psi_{k,\delta}(\alpha) = \frac{\delta}{N}\sum_{m=k}^{N-1} \frac{\binom{m}{k}}{\binom{N}{k}} (1-\alpha)^{-(N-m)} $$
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
    # where I_α is the regularized incomplete beta.
    ratio = T(betainc(k + 1, N - k, αT)) / αT^(k + 1)
    coeff = (one(T) - αT)^(k - N) / T(binomial(N, k))
    return βT / T(N) * coeff * ratio
end

function violation(dist::OneTailViolation, β::Real; tol=1e-10)
    if numsuccess(dist) == 0
        ϵ = 1
    else
        α_lower = 0
        α_upper = 1
        while α_upper - α_lower > tol
            α = (α_lower + α_upper) / 2

            left = β * betainc(numfailure(dist) + 1, numsuccess(dist), α)
            right = α * dist.N * (betainc(numfailure(dist), numsuccess(dist) + 1, α) - betainc(numfailure(dist) + 1, numsuccess(dist), α))
            
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

function confidence(dist::OneTailViolation, ϵ::Tuple{<:Real, <:Real}; tol=1e-10)
    throw(ArgumentError("confidence is not implemented for OneTailViolation"))
end
