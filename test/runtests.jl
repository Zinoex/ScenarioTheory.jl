using ScenarioTheory, StatsFuns, SpecialFunctions
using Test, Supposition, Supposition.Data

@testset "ScenarioTheory.jl" begin
    include("scenario_opt.jl")
    include("wait_and_judge.jl")
    include("onetail.jl")
    include("twotail.jl")
    include("sample_discarding.jl")
    include("enhanced_sample_discarding.jl")
end
