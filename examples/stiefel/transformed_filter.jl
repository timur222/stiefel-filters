include("../../src/stiefel/points.jl")
# import res

τ = 0.8
N = 6
M = 2
eps = .128

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

#println("zero_points = ", zero_points)
#println("infinity_points = ", infinity_points)
#println("coefficient = ", coefficient)
#println("x_max = ", x_max)

# !!! ϵ = 0.12 #??? Параметр исходного преобразования.
x_1_var = x_1(x, y) #??? Нужна для новой нормировки x₁ ↦ -1.
x_2_var = x_2(x, y) #??? Её образ обозначаем α.
x_3_var = x_3(x, y) #??? Первая опорная точка нового преобразования.
x_6_var = x_6(x, y)
x_5_var = x_5(x, y)
x_4_var = x_4(x, y) #??? Вторая опорная точка нового преобразования.

# !!! c = (ϵ/x_4_var - 2/x_6_var)/(ϵ - 2)
# !!! d = 2/x_6_var - c
# !!! c = ((x_4_var + 1) / ϵ - (x_6_var + 1) / 2) / (x_4_var - x_6_var) #??? transform(x₄)=ϵ, transform(x₆)=2.
# !!! c = ((x_4_var + 1) / ϵ - (x_6_var + 1) / 2) / (x_4_var - x_6_var) #??? До сдвига T₀(x₄)=ϵ, T₀(x₆)=2.
# !!! d = (x_4_var + 1) / ϵ - c * x_4_var #??? transform(-1)=0 следует из числителя x+1.
# !!! d = (x_4_var + 1) / ϵ - c * x_4_var #??? До сдвига T₀(-1)=0.
# !!! println("Проверка преобразования: T(-1)=", (-1 + 1) / (-c + d),
# !!!         ", T(x₄)=", (x_4_var + 1) / (c * x_4_var + d),
# !!!         ", T(x₆)=", (x_6_var + 1) / (c * x_6_var + d))



# !!! transform = x -> (x+1)/(c*x+d)
# !!! transform = x -> (x+1)/(c*x+d) - 1 #??? Исходное преобразование.
# !!! Сначала нормируем x₃ ↦ -1 и x₆ ↦ 1, затем используем преобразование,
# !!! фиксирующее ±1: v ↦ (v+k)/(1+k*v). Выбираем ветвь |k|<1.
#= !!! Прежняя нормировка: x₃ ↦ -1, x₆ ↦ 1, T(x₄)=-T(x₅).
reference_points = (x_3_var, x_4_var, x_5_var, x_6_var) #???
println("Используется обновлённый transform из ", @__FILE__) #???
println("Опорные точки transform (x₃, x₄, x₅, x₆) = ", reference_points) #???
all(isfinite, reference_points) && x_3_var != x_6_var || #???
    throw(ArgumentError("для transform нужны конечные x₃, x₄, x₅, x₆ и x₃ ≠ x₆")) #???
#= !!! Прежняя нормировка округляла образы x₄ и x₅ до -1 при очень большом x₆.
mobius_span = x_6_var - x_3_var #???
mobius_normalized = s -> 2 * (s - x_3_var) / mobius_span - 1 #???
v4, v5 = mobius_normalized(x_4_var), mobius_normalized(x_5_var) #???
#= !!! Прежний расчёт в Float64 мог округлить допустимое k почти у ±1 ровно до ±1.
sum_v = v4 + v5 #???
cross_v = 1 + v4 * v5 #???
discriminant = (1 - v4^2) * (1 - v5^2) #???
discriminant >= -64 * Base.eps(max(1.0, abs(cross_v)^2, abs(sum_v)^2)) || #???
    throw(ArgumentError("для этих четырёх точек нет действительного дробно-линейного преобразования с требуемой симметрией")) #???
root_discriminant = sqrt(max(0.0, discriminant)) #???
denominator_k = cross_v + copysign(root_discriminant, cross_v) #???
denominator_k != 0 || throw(ArgumentError("опорные точки задают вырожденное преобразование")) #???
mobius_k = -sum_v / denominator_k #???
abs(mobius_k) < 1 || throw(ArgumentError("преобразование вырождается: |k| должно быть меньше 1")) #???
=#
mobius_k = setprecision(BigFloat, 256) do #??? Подбор одного параметра с запасом точности.
    p3, p4, p5, p6 = BigFloat.(reference_points) #???
    w4 = 2 * (p4 - p3) / (p6 - p3) - 1 #???
    w5 = 2 * (p5 - p3) / (p6 - p3) - 1 #???
    if (w4 == -1 && w5 == 1) || (w4 == 1 && w5 == -1) #???
        return 0.0 #??? Образы уже симметричны; подойдёт аффинная нормировка.
    end #???
    b = 1 + w4 * w5 #???
    s = w4 + w5 #???
    delta = (1 - w4^2) * (1 - w5^2) #???
    delta >= 0 || throw(ArgumentError("для этих точек не существует действительного преобразования с требуемой симметрией: $reference_points")) #???
    denominator = b + copysign(sqrt(delta), b) #???
    denominator != 0 || throw(ArgumentError("вырожденные опорные точки: $reference_points")) #???
    k_big = -s / denominator #??? Корень с |k| ≤ 1 и без полюса между x₃ и x₆.
    abs(k_big) < 1 || throw(ArgumentError("условия несовместимы: x₄ или x₅ совпадает с x₃ или x₆; точки: $reference_points")) #???
    k_float = Float64(k_big) #???
    abs(k_float) < 1 || throw(ArgumentError("опорные точки слишком близки для Float64; нужен более точный расчёт x₃…x₆: $reference_points")) #???
    return k_float #???
end #???
all(v -> 1 + mobius_k * v != 0, (v4, v5)) || #???
    throw(ArgumentError("одна из точек x₄, x₅ попала в полюс преобразования")) #???
transform = s -> (mobius_normalized(s) + mobius_k) / (1 + mobius_k * mobius_normalized(s)) #???
=#
#??? Используем координату w=(s-x₃)/(x₆-s), которая сохраняет различие
#??? между x₄ и x₅ даже при x₆ порядка 10³².
w4 = (x_4_var - x_3_var) / (x_6_var - x_4_var) #???
w5 = (x_5_var - x_3_var) / (x_6_var - x_5_var) #???
# !!! isfinite(w4) && isfinite(w5) && w4 > 0 && w5 > 0 ||
# !!!     throw(ArgumentError("для transform нужны x₄ и x₅ строго между x₃ и x₆"))
#??? Для вещественной симметрии достаточно, чтобы w₄ и w₅ были одного знака.
isfinite(w4) && isfinite(w5) && w4 != 0 && w5 != 0 && signbit(w4) == signbit(w5) || #???
    throw(ArgumentError("для этих опорных точек нет вещественного transform с T(x₄)=-T(x₅): точки=$reference_points, w₄=$w4, w₅=$w5")) #???
# !!! mobius_scale = sqrt(w4) * sqrt(w5)
# !!! mobius_scale = sqrt(abs(w4)) * sqrt(abs(w5)) #??? При w₄,w₅<0 образы выходили за (-1,1).
# !!! mobius_scale = copysign(sqrt(abs(w4)) * sqrt(abs(w5)), w4) #??? Эквивалентная запись выбора ветви через знак.
base_scale = sqrt(abs(w4)) * sqrt(abs(w5)) #??? Сначала строим основное преобразование.
base_transform = s -> ((s - x_3_var) - base_scale * (x_6_var - s)) / #???
                      ((s - x_3_var) + base_scale * (x_6_var - s)) #???
base_image4, base_image5 = base_transform(x_4_var), base_transform(x_5_var) #???
all(v -> !isnan(v), (base_image4, base_image5)) || #???
    throw(ArgumentError("не удалось вычислить образы x₄ и x₅: точки=$reference_points")) #???
invert_transform = abs(base_image4) >= 1 || abs(base_image5) >= 1 #??? Если образ снаружи, заменяем T на 1/T.
transform = s -> invert_transform ? inv(base_transform(s)) : base_transform(s) #???
mobius_scale = invert_transform ? -base_scale : base_scale #??? Тот же выбор ветви для обратной формулы.
isfinite(mobius_scale) && mobius_scale != 0 || #???
    throw(ArgumentError("не удалось вычислить масштаб transform: точки=$reference_points, w₄=$w4, w₅=$w5")) #???
# !!! transform = s -> ((s - x_3_var) - mobius_scale * (x_6_var - s)) /
# !!!                  ((s - x_3_var) + mobius_scale * (x_6_var - s))
transformed_middle_points = (transform(x_4_var), transform(x_5_var)) #???
all(isfinite, transformed_middle_points) || #???
    throw(ArgumentError("x₄ или x₅ попала в полюс transform: точки=$reference_points, w₄=$w4, w₅=$w5")) #???
# !!! all(v -> -1 < v < 1, (transform(x_4_var), transform(x_5_var))) ||
# !!!     throw(ArgumentError("образы x₄ и x₅ не попали строго между -1 и 1: точки=$reference_points, образы=$(transform(x_4_var), transform(x_5_var))"))
all(v -> -1 < v < 1, transformed_middle_points) || #???
    throw(ArgumentError("образы x₄ и x₅ не попали строго между -1 и 1: точки=$reference_points, образы=$transformed_middle_points")) #???
# !!! Формула из slice_mobius_transform не подходит для данного x₆.
# !!! Z = x_6_var
# !!! isfinite(Z) && Z > 0 || throw(ArgumentError("для transform нужно конечное x₆ > 0"))
# !!! root = sqrt(1 + Float64(Z))
# !!! L = Float64(Z) / (root + 1)^2
# !!! a = -Float64(Z) / (root + 1)
# !!! transform = s -> isinf(s) ? inv(L) : (s + a) / (L * s - a)
# !!! transform_inv = x -> (d*x-1)/(c-x)
# !!! transform_inv = u -> (d*u-1)/(1-c*u)
# !!! transform_inv = u -> (d*(u+1)-1)/(1-c*(u+1)) #??? Обратное к исходному transform.
# !!! transform_inv = u -> x_3_var + mobius_span * ((u - mobius_k) / (1 - mobius_k * u) + 1) / 2
transform_inv = u -> ((1 - u) * x_3_var + mobius_scale * (1 + u) * x_6_var) / #???
                       ((1 - u) + mobius_scale * (1 + u)) #???
# !!! transform_inv = u -> a * (u + 1) / (L * u - 1)
=#
normalization = symmetric_mobius_normalization(x_1_var, x_4_var, x_5_var, x_6_var) #???
transform = normalization.transform #??? x₁ ↦ -1, x₄ ↦ 1, x₅ ↦ 1/γ, x₆ ↦ -1/γ.
transform_inv = normalization.inverse #??? Обратное к этой же нормировке.
α = transform(x_2_var) #???
β = transform(x_3_var) #???
γ = normalization.gamma #???
normalization_path = joinpath(@__DIR__, "../../data/normalization_points.txt") #??? Те же числа использует сечение.
open(normalization_path, "w") do io #???
    for point in (x_1_var, x_2_var, x_3_var, x_4_var, x_5_var, x_6_var) #???
        println(io, repr(Float64(point))) #???
    end #???
end #???
println("Опорные точки сохранены для сечения: ", normalization_path) #???
println("Нормировка: T(x₁)=", transform(x_1_var), ", T(x₂)=α=", α, #???
        ", T(x₃)=β=", β, ", T(x₄)=", transform(x_4_var), #???
        ", T(x₅)=1/γ=", transform(x_5_var), ", T(x₆)=-1/γ=", transform(x_6_var), #???
        ", γ=", γ) #???
# !!! println("Проверка преобразования: T(-1)=", transform(-1),
# !!!         ", T(x₄)=", transform(x_4_var), ", T(x₆)=", transform(x_6_var))
# !!! println("Проверка преобразования: T(x₃)=", transform(x_3_var),
# !!!         ", T(x₆)=", transform(x_6_var),
# !!!         ", T(x₄)=", transform(x_4_var), ", T(x₅)=", transform(x_5_var),
# !!!         ", T(x₄)+T(x₅)=", transform(x_4_var) + transform(x_5_var))

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

#??? Проверяем именно значения кривой, а не горизонтальных полос.
# !!! for (name, q) in (("-1", -1.0), ("0", 0.0), ("x₄", x_4_var),
# !!! for (name, q) in (("x₃", x_3_var), ("0", 0.0), ("x₄", x_4_var),
# !!!                   ("x₅", x_5_var), ("x₆", x_6_var))
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

println("transform_inv(x_4_var) / transform_inv(x_5_var) = ", transform_inv(x_4_var)/transform_inv(x_5_var))

println("OM[1, 1] = ", OM[1, 1])
println("OM[2, 2] = ", OM[2, 2])

# !!! x_1_var = x_1(x, y) #??? Уже вычислено перед transform.
# !!! x_2_var = x_2(x, y) #??? Уже вычислено перед transform.
# !!! x_3_var = x_3(x, y) #??? Вычисляется выше, до определения transform.
# !!! x_4_var = x_4(x, y) #??? Уже вычислено перед transform.
# !!! x_5_var = x_5(x, y) #??? Уже вычислено перед transform.
# !!! x_6_var = x_6(x, y) #??? Уже вычислено перед transform.

x_4_prime_var = x_4_prime(x, y)
x_5_prime_var = x_5_prime(x, y)

graph_name = "up_4" #??? Измените эту строку, чтобы назвать график и PNG-файл.
# !!! draw_graph(frac_real, -4, 4;
plot_radius = max(2.0, 1.1 * maximum(abs, transform.((x_1_var, x_2_var, x_3_var, #???
                                                       x_4_var, x_5_var, x_6_var)))) #???
# !!! fig = draw_graph(frac_real, -2, 2100; #??? Диапазон предыдущей нормировки.
fig = draw_graph(frac_real, -plot_radius, plot_radius; #???
    n=1_000_000,
    ylims=(-1.5, 1.5), #??? Пределы оси y; здесь можно задать нужный диапазон.
    # !!! intervals=((real(-1000_000), real(transform(-1))),
    # !!! (real(transform(0)), real(transform(x_4_var))),
    # !!! (real(transform(x_5_var)), real(transform(x_6_var))),
    # !!! (real(transform(-1_000_000)), real(100))),
    # !!! intervals=((real(-1_000_000), real(transform(x_2(x, y)))),
    # !!!            (real(-1), real(transform(x_4_var))),
    # !!!            (real(transform(x_5_var)), real(1)),
    # !!!            (real(transform(x_1(x, y))), real(100))),
    intervals=((-1.0, α), #??? x₁, x₂ → -1, α.
               (β, 1.0), #??? x₃, x₄ → β, 1.
               (inv(γ), -inv(γ)), #??? x₅, x₆ → 1/γ, -1/γ.
               (1.0, Inf)), #??? Четвёртая красная полоса по выбору пользователя.
    
    values=(-(1.0+result.value)/2, (1.0+result.value)/2, (1.0+result.value)/2, -(1.0+result.value)/2))

# !!! Заголовок и имя файла берутся из одного настраиваемого названия.
# !!! fig.content[1].title = graph_name
fig.content[1].title = "" #??? Перед сохранением явно убираем верхний заголовок.
# !!! graph_path = joinpath(@__DIR__, "test100_graph.png")
graph_path = joinpath(mkpath(joinpath(@__DIR__, "../../output")), graph_name * ".png") #???
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
