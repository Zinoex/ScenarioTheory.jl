using CairoMakie
using ScenarioTheory

N = 10000
beta = 1e-6

k_values = 0:10:N


eps_values = [last(violation(OneTailViolation(k, N), beta)) for k in k_values]
epsLU_values = [violation(TwoTailViolation(k, N), beta) for k in k_values]

fig = Figure()
ax = Axis(fig[1, 1], xlabel="k", ylabel="ε", title="ε vs k")
lines!(ax, k_values, eps_values, label="ε")
lines!(ax, k_values, [epsLU[1] for epsLU in epsLU_values], label="ε_L")
lines!(ax, k_values, [epsLU[2] for epsLU in epsLU_values], label="ε_U")
axislegend(ax)

display(fig)