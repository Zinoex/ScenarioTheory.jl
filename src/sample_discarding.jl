struct SampleDiscarding <: AbstractScenarioProblem
    samples::Int
    decision_vars::Int
    discarded::Int

    function SampleDiscarding(samples::Int, decision_vars::Int, discarded::Int)
         if samples < 1
            throw(DomainError(samples, "expected samples ≥ 1"))
        end

        if discarded < 0
            throw(DomainError(discarded, "expected discarded ≥ 0"))
        end

        if decision_vars < 1
            throw(DomainError(decision_vars, "expected decision_vars ≥ 1"))
        end
        
        if decision_vars + discarded > samples
            throw(ArgumentError("expected decision_vars + discarded ≤ samples"))
        end
        
        return new(samples, decision_vars, discarded)
    end
end

function violation(dist::SampleDiscarding, β::Real; tol=1e-10)
    N = dist.samples
    d = dist.decision_vars
    k = dist.discarded
    r = k + d - 1

    if N == d
        ϵ = 1.0
    else
        α_lower = 0.0
        α_upper = 1.0

        while α_upper - α_lower > tol
            α = (α_lower + α_upper) / 2
            β_mid = binomial(r, k) * binomcdf(N, r, α)

            if β_mid > β
                α_upper = α
            else
                α_lower = α
            end
        end

        ϵ = α_upper
    end

    return zero(ϵ), ϵ
end
