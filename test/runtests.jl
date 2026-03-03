using ScenarioTheory, StatsFuns
using Test, Supposition, Supposition.Data

@testset "ScenarioTheory.jl" begin
    include("scenario_opt.jl")
    include("wait_and_judge.jl")
    include("onetail.jl")
    include("twotail.jl")
end
