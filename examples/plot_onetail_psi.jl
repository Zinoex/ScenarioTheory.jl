using CairoMakie
using StatsFuns

function eq(β, N, k)
    function fun(α)
        return β * binomccdf(N, α, k) - α * N * binompdf(N, α, k)
    end

    return fun
end

β = 1e-6
N = 1000

fig = Figure()
ax = Axis(fig[1, 1], xlabel="α", ylabel="ψ(α)", title="One-tail ψ(α) for β=$β, N=$N")

k = 200
f = eq(β, N, k)
lines!(ax, 0:0.001:1, f, label="k=$k")

k = 210
f = eq(β, N, k)
lines!(ax, 0:0.001:1, f, label="k=$k")

# hlines!(ax, [0], color=:red, linestyle=:dash, label="ψ(α)=1")
axislegend(ax, position=:rb)

display(fig)