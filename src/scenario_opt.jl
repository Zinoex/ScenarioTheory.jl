struct ScenarioOptimization <: AbstractScenarioTheory
    samples::Int
    decision_vars::Int

    function ScenarioOptimization(samples::Int, decision_vars::Int)
        if samples < 1
            throw(DomainError(samples, "expected samples ≥ 1"))
        end

        if decision_vars < 1
            throw(DomainError(decision_vars, "expected decision_vars ≥ 1"))
        end
        
        if decision_vars > samples
            throw(ArgumentError("expected decision_vars ≤ samples"))
        end
        
        return new(samples, decision_vars)
    end
end

function violation(dist::ScenarioOptimization, β::Real; α_tol=1e-10, float_tol=20)
    N = dist.samples
    d = dist.decision_vars
    k = d - 1

    if N == d
        ϵ = 1.0
    else
        α_lower = 0.0
        α_upper = 1.0

        while α_upper > nextfloat(α_lower)
            α = (α_lower + α_upper) / 2
            β_mid = binomcdf(N, α, k)

            if β > nextfloat(β_mid, float_tol)
                α_upper = α
            else
                α_lower = α
            end
        end

        ϵ = α_upper
    end

    return zero(ϵ), ϵ
end
