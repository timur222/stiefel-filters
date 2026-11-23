# !!! include("../../src/stiefel/func.jl")
include("../../src/stiefel/points.jl")
# !!! include("../../src/stiefel/points.jl")
# import res

τ = 0.2
N = 6
# !!! p = [im * τ]
p = ComplexF64[im * τ * N ;;] #??? Матрица периода рода 1 должна иметь размер 1×1.

zero_points, infinity_points = shtiff_fraction_gen_1(τ, N, p; tol=1e-10)

frac = fraction(zero_points, infinity_points, -1)
frac_real = x -> real(frac(x))

#transform(x)

println("here")
draw_graph(frac_real, -10, 10; n=1_000_000)

println("x_1 ", x_gen_1(τ, N; special_point = 1))
println("x_2 ", x_gen_1(τ, N; special_point = 2))
println("x_3 ", x_gen_1(τ, N; special_point = 3))
println("x_4 ", x_gen_1(τ, N; special_point = 4))

println("x_max 1 ", x_gen_1(τ, N; right_point = 1))
println("x_max 2 ", x_gen_1(τ, N; right_point = 2))
println("x_max 3 ", x_gen_1(τ, N; right_point = 3))
