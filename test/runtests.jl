using ScenarioTheory
using Test, Supposition, Supposition.Data

@testset "ScenarioTheory.jl" begin
    include("scenario_opt.jl")
    include("wait_and_judge.jl")
end
