module ScenarioTheory

using SpecialFunctions: beta_inc
betainc(a, b, x) = first(beta_inc(a, b, x))

abstract type AbstractViolation end
numsuccess(dist::AbstractViolation) = numscenarios(dist) - numfailure(dist)

export psi, violation, confidence

include("onetail.jl")
export OneTailViolation

include("twotail.jl")
export TwoTailViolation

end
