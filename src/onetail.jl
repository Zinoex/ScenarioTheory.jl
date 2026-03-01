struct CompresionOneTail <: AbstractScenarioProblem
    samples::Int
    compressed::Int

    function CompresionOneTail(samples::Int, compressed::Int)
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

function psi(dist::CompresionOneTail, β::Real, α::Real)
    # $$ \Psi_{k,\delta}(\alpha) = \frac{\delta}{N}\sum_{m=k}^{N-1} \frac{\binom{m}{k}}{\binom{N}{k}} (1-\alpha)^{-(N-m)} $$
    k = dist.compressed
    N = dist.samples

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

function violation(dist::CompresionOneTail, β::Real; tol=1e-10)
    N = dist.samples
    k = dist.compressed
    l = N - k

    if l == 0
        ϵ = 1.0
    else
        α_lower = 0.0
        α_upper = 1.0
        while α_upper - α_lower > tol
            α = (α_lower + α_upper) / 2

            left = β * betainc(k + 1, l, α)
            right = α * N * (betainc(k, l + 1, α) - betainc(k + 1, l, α))
            
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
