module ScenarioTheory

using StatsFuns: binompdf, binomcdf, binomccdf

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
