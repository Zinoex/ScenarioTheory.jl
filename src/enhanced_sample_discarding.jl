struct EnhancedSampleDiscarding <: AbstractScenarioTheory
    samples::Int
    decision_vars::Int
    discarding_rounds::Int

    function EnhancedSampleDiscarding(samples::Int, decision_vars::Int, discarding_rounds::Int)
         if samples < 1
            throw(DomainError(samples, "expected samples ≥ 1"))
        end

        if decision_vars < 1
            throw(DomainError(decision_vars, "expected decision_vars ≥ 1"))
        end

        if discarding_rounds < 0
            throw(DomainError(discarding_rounds, "expected discarding_rounds ≥ 0"))
        end
        
        if decision_vars * (1 + discarding_rounds) > samples
            throw(ArgumentError("expected decision_vars * (1 + discarding_rounds) ≤ samples"))
        end
        
        return new(samples, decision_vars, discarding_rounds)
    end
end

function violation(dist::EnhancedSampleDiscarding, β::Real; tol_steps=20)
    N = dist.samples
    d = dist.decision_vars
    k = dist.discarding_rounds
    r = k + d - 1

    if N == d * (1 + k)
        ϵ = 1.0
    else
        α_lower = 0.0
        α_upper = 1.0

        while α_upper > nextfloat(α_lower)
            α = (α_lower + α_upper) / 2
            β_mid = binomcdf(N, α, r)

            if β > nextfloat(β_mid, tol_steps)
                α_upper = α
            else
                α_lower = α
            end
        end

        ϵ = α_upper
    end

    return zero(ϵ), ϵ
end
