module ScenarioTheory

using SpecialFunctions: beta_inc
betainc(a, b, x) = first(beta_inc(a, b, x))

abstract type AbstractScenarioProblem end

export violation

include("onetail.jl")
export CompressionOneTail

include("twotail.jl")
export CompressionTwoTail

include("scenario_opt.jl")
export ScenarioOptimization

include("wait_and_judge.jl")
export WaitAndJudgeScenarioOptimization

include("sample_discarding.jl")
export SampleDiscardingScenarioOptimization

end
