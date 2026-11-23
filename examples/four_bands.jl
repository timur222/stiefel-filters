include("../src/four_bands/StiefelFourBands.jl")
using .StiefelFourBands

# Высота клетки, степень дроби и два уровня разрезов, считая снизу.
τ = 0.7000585045606659
N = 3
M1, M2 = 1, 2

# В используемом базисе Ω[1,:] = τ/2 * [N, N-M1, N-M2].
# Две диагональные компоненты свободны, Ω₂₃ получается из кривой.
filter_model = filter_from_periods(τ, N, M1, M2;
    omega22=0.8522154941218228, omega33=0.5744824994896586,
    trims=(0.2005990182541181, 0.30045953696537414))

p = im * period_parameters(filter_model).omega
println("Матрица периодов p = ", p)
println("Концы четырёх полос = ", filter_model.endpoints)
curve = ThetaCurve(filter_model.geom)
frac(x) = evaluate(filter_model, x; method=:theta, curve).direct
println("f(0.3 + 0.2i) = ", frac(0.3 + 0.2im))

output_dir = mkpath(joinpath(@__DIR__, "../output"))
println("График: ", plot_filter(joinpath(output_dir, "four_bands.html"), filter_model))
