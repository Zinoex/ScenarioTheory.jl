struct WaitAndJudgeScenarioOptimization <: AbstractScenarioProblem
    samples::Int
    support::Int

    function WaitAndJudgeScenarioOptimization(samples::Int, support::Int)
        if samples < 1
            throw(DomainError(samples, "expected samples ≥ 1"))
        end

        if support < 1
            throw(DomainError(support, "expected support ≥ 1"))
        end
        
        if support > samples
            throw(ArgumentError("expected support ≤ samples"))
        end
        
        return new(samples, support)
    end
end

function violation(dist::WaitAndJudgeScenarioOptimization, β::Real; tol=1e-10)
    N = dist.samples
    k = dist.support
    l = N - k

    if N == k
        ϵ = 1.0
    else
        α_lower = 0.0
        α_upper = 1.0
        while α_upper - α_lower > tol
            α = (α_lower + α_upper) / 2

            # \frac{\beta}{N + 1} \sum_{m=k}^N \binom{m}{k}(1-\alpha)^{m - k} - \binom{N}{k}(1-\alpha)^{N - k} = 0

            # \frac{\beta}{N + 1} \sum_{m=k}^N \binom{m}{k}(1-\alpha)^{m - k} = \binom{N}{k}(1-\alpha)^{N - k}

            # \frac{\beta}{N + 1} \frac{(1 - \sum_{i=0}^k \binom{N + 1}{i}\alpha^{i}(1-\alpha)^{N + 1 - i}}{\alpha^{k + 1}} = \binom{N}{k}(1-\alpha)^{N - k}

            # \frac{\beta}{N + 1} (1 - \sum_{i=0}^k \binom{N + 1}{i}\alpha^{i}(1-\alpha)^{N + 1 - i} = \binom{N}{k}\alpha^{k + 1}(1-\alpha)^{N - k}

            # \beta (1 - \sum_{i=0}^k \binom{N + 1}{i}\alpha^{i}(1-\alpha)^{N + 1 - i} = \alpha * (N + 1) * \binom{N}{k}\alpha^k(1-\alpha)^{N - k}

            left = β * binomccdf(N + 1, α, k)
            right = α * (N + 1) * binompdf(N, α, k)
            
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
