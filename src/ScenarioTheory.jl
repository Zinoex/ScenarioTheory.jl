module ScenarioTheory

using StatsFuns: binompdf, binomcdf, binomccdf, binomlogcdf
using SpecialFunctions: logabsbinomial

abstract type AbstractScenarioTheory end

export violation

include("onetail.jl")
export CompressionOneTail

include("twotail.jl")
export CompressionTwoTail

include("scenario_opt.jl")
export ScenarioOptimization

include("wait_and_judge.jl")
export WaitAndJudge

include("sample_discarding.jl")
export SampleDiscarding

include("enhanced_sample_discarding.jl")
export EnhancedSampleDiscarding

end
