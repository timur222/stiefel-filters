using CairoMakie
using Printf
include("../../src/four_bands/StiefelFourBands.jl")
using .StiefelFourBands
include("samples.jl")

# julia --project=. examples/three_bands/save_function.jl 16 0.1 1 1 0.5
n = isempty(ARGS) ? 16 : parse(Int, ARGS[1])
epsilon = length(ARGS) < 2 ? 0.1 : parse(Float64, ARGS[2])
weights = length(ARGS) <= 2 ? (1.0, 1.0, 0.5) : parse.(Float64, ARGS[3:end])
sol = oscillating_threebands(n, epsilon; weights)
rows = [sample_threeband(sol, band) for band in 1:3]
colors = ("#254d76", "#167459", "#a34526")
labels = ("E₁", "E₂", "E₃")
number(x) = @sprintf("%.4g", x)

fig = Figure(size=(1160, 640), fontsize=16, figure_padding=22)
Label(fig[0, 1:3], "Рациональная функция степени $n", fontsize=25, halign=:left)
Label(fig[1, 1:3], "Веса  ($(join(number.(weights), ", ")))     ·     ε = $(number(epsilon))",
      color=:gray35, halign=:left)
overview = Axis(fig[2, 1:3], xlabel="x", ylabel="f(x)", xticks=-1:0.25:1)
for band in 1:3
    data = rows[band]
    lines!(overview, Float64.(data.x), data.value; color=colors[band], linewidth=2)
    center = band == 1 ? -1 : 1
    width = epsilon / weights[band]
    levels = [center-width, center, center+width]
    a, b = sol.ends[2band-1:2band]
    ax = Axis(fig[3, band], title="$(labels[band])   [$(number(a)), $(number(b))]",
              titlesize=15, xlabel="log₁₀|x|", ylabel=band==1 ? "f(x)" : "",
              yticks=(levels, number.(levels)))
    hlines!(ax, levels[[1,3]]; color=:gray65, linestyle=:dash, linewidth=1)
    lines!(ax, data.logx, data.value; color=colors[band], linewidth=2)
    xlims!(ax, first(data.logx), last(data.logx))
    ylims!(ax, center-1.18width, center+1.18width)
end
xlims!(overview, -1, 1)
rowgap!(fig.layout, 8)
out = mkpath(joinpath(@__DIR__, "../../output"))
for extension in ("svg", "png")
    path = joinpath(out, "rational-function-n$n.$extension")
    save(path, fig)
    println(abspath(path))
end
