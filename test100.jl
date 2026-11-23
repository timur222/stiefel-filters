include("func.jl")
include("points.jl")
# import res

τ = 0.7
N = 10
M = 4
eps = 0.06
x1 = -2
x2 = 2
graph_name = "down_3"

# размеры
H1 = N * τ
H2 = 2
H3 = (N - M) * τ
H4 = 0
H5 = M * τ

Ω₁₁ = H1 / H2
Ω₂₁ = H5 / H2
Ω₁₂ = H5 / H2
Ω₂₂ = H5 / H2

OM = [Ω₁₁ Ω₁₂; Ω₂₁ Ω₂₂]
OM[2, 2] *= 1 + eps

x = OM[1, 1]
y = OM[2, 2] - x * M / N  


println(OM)
p = OM * im

t = τ
n = N
m = M

e = special_chars()
u = special_points(p)

points_m = (m * t * im - 1) / (-4)
res_new = shtiff_fraction_new(t, n, p, m)
zero_points, infinity_points, coefficient, x_max = res_new[1], res_new[2], res_new[3], 4 * res_new[4]

x_1_var = x_1(x, y) #??? Нужна для новой нормировки x₁ ↦ -1.
x_2_var = x_2(x, y) #??? Её образ обозначаем α.
x_3_var = x_3(x, y) #??? Первая опорная точка нового преобразования.
x_6_var = x_6(x, y)
x_5_var = x_5(x, y)
x_4_var = x_4(x, y) #??? Вторая опорная точка нового преобразования.



# !!! normalization = symmetric_mobius_normalization(x_1_var, x_4_var, x_5_var, x_6_var)
normalization = six_point_mobius_normalization(x_1_var, x_2_var, x_3_var, #???
                                                x_4_var, x_5_var, x_6_var) #???
# !!! transform = normalization.transform #??? x₁ ↦ -1, x₄ ↦ 1, x₅ ↦ 1/γ, x₆ ↦ -1/γ.
transform = normalization.transform #??? x₁,x₂ ↦ ±1/γ; x₃ ↦ 1; x₆ ↦ -1.
transform_inv = normalization.inverse #??? Обратное к этой же нормировке.
# !!! α = transform(x_2_var)
# !!! β = transform(x_3_var)
α = normalization.alpha #??? Образ x₅.
β = normalization.beta #??? Образ x₄.
γ = normalization.gamma #???
first_interval_wraps = min(x_1_var, x_2_var) < normalization.pole < #???
                       max(x_1_var, x_2_var) #???
all(v -> isfinite(v) && abs(v) < 1, (α, β)) || #???
    throw(ArgumentError("нужны |α|<1 и |β|<1; получено α=$α, β=$β для x₁…x₄=$((x_1_var, x_2_var, x_3_var, x_4_var))")) #???
normalization_path = joinpath(@__DIR__, "test100_normalization_points.txt") #??? Те же числа использует сечение.
open(normalization_path, "w") do io #???
    for point in (x_1_var, x_2_var, x_3_var, x_4_var, x_5_var, x_6_var) #???
        println(io, repr(Float64(point))) #???
    end #???
end #???
println("Опорные точки сохранены для сечения: ", normalization_path) #???
# !!! println("Нормировка: T(x₁)=", transform(x_1_var), ", T(x₂)=α=", α,
# !!!         ", T(x₃)=β=", β, ", T(x₄)=", transform(x_4_var),
# !!!         ", T(x₅)=1/γ=", transform(x_5_var), ", T(x₆)=-1/γ=", transform(x_6_var), ", γ=", γ)
println("Нормировка: T(x₁)=1/γ=", transform(x_1_var), #???
        ", T(x₂)=-1/γ=", transform(x_2_var), ", T(x₃)=", transform(x_3_var), #???
        ", T(x₄)=β=", β, ", T(x₅)=α=", α, #???
        ", T(x₆)=", transform(x_6_var), ", γ=", γ) #???

frac = fraction(zero_points, infinity_points, coefficient)
frac_1 = x -> real(frac(x))
frac_2 = x -> -1*real(frac(x))
result = find_maximum(frac_1, 0, 4; n=100_001)
result2 = find_maximum(frac_2, 0, 4; n=100_001)

#??? Уточняем граничные точки по самой рациональной функции: Re(frac(q)) = ±1.
function level_boundary_near(frac, q, left, right; target=nothing) #???
    original_q = q #???
    value = real(frac(q)) #???
    isfinite(value) || throw(ArgumentError("значение дроби в исходной границе $q не конечно")) #???
    level = isnothing(target) ? (abs(value - 1) <= abs(value + 1) ? 1.0 : -1.0) : target #???
    residual = abs(value - level) #???
    residual <= 1e-10 && return q #??? Точная граничная точка не требует сдвига.
    score(s) = begin #???
        v = real(frac(s)) #???
        isfinite(v) ? -abs(v - level) : -1e300 #??? Полюс не может быть границей уровня ±1.
    end #???
    candidate = find_maximum(score, left, right; n=2_001, iterations=48).x #???
    candidate_residual = abs(real(frac(candidate)) - level) #???
    edge_margin = 32 * Base.eps(max(abs(left), abs(right), 1.0)) #???
    if isfinite(candidate_residual) && candidate_residual <= 1e-5 && #???
       candidate_residual < residual && left + edge_margin < candidate < right - edge_margin #???
        q, residual = candidate, candidate_residual #???
    end #???
    residual > 1e-5 && @warn "Не удалось уточнить границу до уровня ±1; оставлена исходная точка" original_q q level residual left right #???
    return q #???
end #???

# !!! nominal_boundaries = (-1.0, 0.0, x_4_var, x_5_var, x_6_var)
nominal_boundaries = (x_1_var, x_2_var, x_4_var, x_5_var, x_6_var) #??? Границы новой нормировки.
all(isfinite, nominal_boundaries) && length(unique(nominal_boundaries)) == 5 || #???
    throw(ArgumentError("граничные точки должны быть конечными и различными")) #???
refined_boundaries = ntuple(i -> begin #???
    q = nominal_boundaries[i] #???
    lower = filter(v -> v < q, nominal_boundaries) #???
    upper = filter(v -> v > q, nominal_boundaries) #???
    left = isempty(lower) ? q - (minimum(upper) - q) / 2 : (maximum(lower) + q) / 2 #???
    right = isempty(upper) ? q + (q - maximum(lower)) / 2 : (q + minimum(upper)) / 2 #???
    level_boundary_near(frac, q, left, right; target=i == 1 ? -1.0 : (i == 2 ? 1.0 : nothing)) #???
end, 5) #???
# !!! u_minus_one, u_zero, u_x4, u_x5, u_x6 = transform.(refined_boundaries)
u_x1, u_x2, u_x4, u_x5, u_x6 = transform.(refined_boundaries) #???
if u_x5 == u_x6 #??? Уточнение не должно уничтожать третий участок.
    original_u_x5, original_u_x6 = transform(x_5_var), transform(x_6_var) #???
    if !isnan(real(original_u_x5)) && !isnan(real(original_u_x6)) && original_u_x5 != original_u_x6 #???
        @warn "Уточнение свело границы третьей полосы; возвращены исходные образы x₅ и x₆" #???
        u_x5, u_x6 = original_u_x5, original_u_x6 #???
    end #???
end #???
println("Уточнённые границы полос: ", refined_boundaries) #???
println("Проверка уровней на границах: ", real.(frac.(refined_boundaries))) #???
println("Границы полос после transform: ", (u_x1, u_x2, u_x4, u_x5, u_x6)) #???

# !!! frac_real = u -> real(frac(u)) #??? Не учитывало смену горизонтальной координаты.
frac_real = u -> real(frac(transform_inv(u))) #??? График исходной дроби в новой координате.

for (name, q) in (("x₁", x_1_var), ("x₂", x_2_var), ("x₃", x_3_var), #???
                  ("x₄", x_4_var), ("x₅", x_5_var), ("x₆", x_6_var)) #???
    u = real(transform(q)) #???
    q_back = transform_inv(u) #???
    original = real(frac(q)) #???
    transformed = frac_real(u) #???
    println("Граница ", name, ": transform(q) = ", u, #???
            ", frac(q) = ", original, ", frac_real(transform(q)) = ", transformed, #???
            ", |transform_inv(transform(q))-q| = ", abs(q_back-q), #???
            ", расстояние до ±1 = ", min(abs(transformed-1), abs(transformed+1))) #???
end #???

#transform(x)

println("OM[1, 1] = ", OM[1, 1])
println("OM[2, 2] = ", OM[2, 2])

x_4_prime_var = x_4_prime(x, y)
x_5_prime_var = x_5_prime(x, y)

# !!! fig = draw_graph(frac_real, -2, 2100; #??? Диапазон предыдущей нормировки.
fig = draw_graph(frac_real, x1, x2; #???
    n=1_000_000,
    ylims=(-2, 2), 
    # !!! intervals=((-1.0, α), (β, 1.0), (inv(γ), -inv(γ)), (1.0, Inf)),
    intervals=((-1/γ, 1/γ), #??? Три полосы ±1/γ.
               (transform(x_3(x, y)), transform(x_4(x, y))), #??? x₃, x₄ → 1, β.
               (transform(x_5(x, y)), transform(x_6(x, y))),
               (-10, -14)), #??? Четвёртая красная полоса по выбору пользователя.
    # !!! wrapped_intervals=(3,), #??? Прежний третий участок проходил через полюс.
    wrapped_intervals=first_interval_wraps ? (1,) : (), #??? Первый участок разделяем, только если он проходит через полюс.
    values=(-(1.0+result.value)/2, (1.0+result.value)/2, (1.0+result.value)/2, -(1.0+result.value)/2))

# !!! Заголовок и имя файла берутся из одного настраиваемого названия.
# !!! fig.content[1].title = graph_name
fig.content[1].title = "" #??? Перед сохранением явно убираем верхний заголовок.
# !!! graph_path = joinpath(@__DIR__, "test100_graph.png")
graph_path = joinpath(@__DIR__, graph_name * ".png") #???
println("Сохраняю рисунок без верхнего заголовка: ", graph_path) #???
Makie.save(graph_path, fig) #???
println("График сохранён: ", graph_path) #???

# !!! Маркеры показывают фактическое значение дроби точно на границах интервалов.
# !!! boundary_q = (-1.0, 0.0, x_4_var, x_5_var, x_6_var)
# !!! boundary_u = [real(transform(q)) for q in boundary_q]
# !!! boundary_f = frac_real.(boundary_u)
# !!! Makie.scatter!(fig.content[1], boundary_u, boundary_f; color=:black, markersize=12)
# !!! fig.content[1].title = "Чёрные точки — значения дроби на границах"

println(transform(-100000))
println(transform(-1))
println(transform(0))
println(transform(x_4_var))
println(transform(x_5_var))
println(transform(x_6_var))
println()
println(transform_inv(-100000))
println(transform_inv(-1))
println(transform_inv(0))
println(transform_inv(x_4_var))
println(transform_inv(x_5_var))
println(transform_inv(x_6_var))

println()
println(transform_inv(-100000))
println(transform_inv(-1))
println(transform_inv(0))
println(transform_inv(x_4_prime(x, y)))
println(transform_inv(x_5_prime(x, y)))
println(transform_inv(x_6(x, y)))


#println("transform(-1) = ", transform(-1))
#println("transform(x_4_var) = ", transform(x_4_var))
#println("transform(x_5_var) = ", transform(x_5_var))
#println("transform(x_6_var) = ", transform(x_6_var))

println("transform_inv(x_4_var) / transform_inv(x_5_var) = ", transform_inv(x_4_var)/transform_inv(x_5_var))

println("OM[1, 1] = ", OM[1, 1])
println("OM[2, 2] = ", OM[2, 2])


println("x_1 ", x_1(x, y))
println("x_2 ", x_2(x, y))
println("x_3 ", x_3(x, y))
println("x_4 ", x_4(x, y))
println("x_4_prime ", x_4_prime(x, y))
println("x_5 ", x_5(x, y))
println("x_5_prime ", x_5_prime(x, y))
println("x_6 ", x_6(x, y))

println("real(transform(transform(x_1_var))) = ", real(transform(x_1(x, y))))
println("real(transform(transform(x_2_var))) = ", real(transform(x_2(x, y))))

println("x_2     = ", x_2(x, y))
println("x_2_new = ", x_2_new(x, y))

println("x_3     = ", x_3(x, y))
println("x_3_new = ", x_3_new(x, y))

println("x_4     = ", x_4(x, y))
println("x_4_new = ", x_4_new(x, y))

println("x_5     = ", x_5(x, y))
println("x_5_new = ", x_5_new(x, y))

println("x_6     = ", x_6(x, y))
println("x_6_new = ", x_6_new(x, y))
