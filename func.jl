using LinearAlgebra
using Polynomials
import Makie
println("Loading theta.jl ...")
flush(stdout)
include("theta.jl")
println("theta.jl loaded.")
flush(stdout)

newton_number = 1000
use_faster = USE_FASTER_THETA

function special_points(p)
    # по основным тэта-характеристикам получить соответствующую точку в якобиане
    a = special_chars()
    b = []
    for i in 1:6
        push!(b, th_char(a[i][:, 1], a[i][:, 2], p))
    end
    return b
end

function th_char(e1, e2, p)
    # столбцы тэта-характеристики -> точка якобиана
    return 0.5 * (p * e1 + e2)
end

function th_char_new(e, p)
    # матрица тэта-характеристики -> точка на якобиане
    e1, e2 = e_columns(e)
    return th_char(e1, e2, p)
end

function th_char_from_u(u, p)
    # точка на якобиане -> тэта-характеристика
    e1 = 2 * real.(u)
    #println("e2 ", e2)
    e2 = 2 * inv(imag(p)) * (imag(u))
    #println("e2 ", e2)
    return vcat(e2', e1')'
end

function eps1()
    # машинный ноль
    eps_c = 1
    while 1 + eps_c > 1
        eps_c /= 2
    end
    return eps_c
end


step_num = 100
newton_steps_count = 5
starting_point_index = 6

function newton(f, fp, x0; tol=1e-10, maxiter=100)
    x = x0

    iter = 0
    diff = Inf
    diff_abs = Inf

    while (iter < maxiter) && (diff > tol)
        iter += 1
        x_new = x - f(x) / fp(x)
        diff_abs = max(abs(x_new - x) / abs(x_new), abs(f(x_new)))
        diff = abs(x_new - x)
        x = x_new
    end

    #println("diff = ", diff)
    #println("diff_abs = ", diff_abs)

    return x, diff_abs, iter, diff
end

function newton(f, jac, x0::AbstractVector; tol=1e-10, maxiter=100)
    x = copy(x0)

    iter = 0
    diff_abs = Inf
    diff = Inf
    diff_val = Inf

    fx = f(x)
    diff_val = norm(fx)
    if diff_val <= tol
        return x, diff_abs, iter, diff, diff_val
    end

    while (iter < maxiter) && (diff > tol) && (diff_val > tol)
        iter += 1

        jx = jac(x)
        delta = jx \ fx
        x_new = x - delta

        fx = f(x_new)
        diff_abs = max(norm(delta) / max(norm(x_new), 1.0), norm(fx))
        diff = norm(delta)
        diff_val = norm(fx)
        x = x_new
    end

    #println("diff = ", diff)
    #println("diff_abs = ", diff_abs)


    return x, diff_abs, iter, diff, diff_val
end

function project_to_image(u, p; shtiff_chars=[[]], act_vars=[1], act_init=nothing, tol=1e-10, maxiter=100)
    s = length(u)
    act_vars = collect(act_vars)
    shtiff_chars = collect(shtiff_chars)

    if isempty(shtiff_chars) || all(isempty, shtiff_chars)
        return copy(u)
    end

    if isempty(act_vars)
        error("project_to_image: act_vars не может быть пустым, если есть уравнения.")
    end

    if length(unique(act_vars)) != length(act_vars)
        error("project_to_image: в act_vars есть повторяющиеся индексы.")
    end

    if any(j -> (j < 1 || j > s), act_vars)
        error("project_to_image: индексы в act_vars должны лежать в 1:length(u).")
    end

    if length(act_vars) != length(shtiff_chars)
        error("project_to_image: для метода Ньютона число активных переменных должно совпадать с числом условий theta = 0.")
    end

    u_start = copy(u)

    function build_u(v)
        u_curr = copy(u_start)
        for (local_index, global_index) in enumerate(act_vars)
            u_curr[global_index] = v[local_index]
        end
        return u_curr
    end

    function f_proj(v)
        u_curr = build_u(v)
        return ComplexF64[
            θ_shtiff(u_curr, p; char_shtiff=char_k, faster=use_faster)
            for char_k in shtiff_chars
        ]
    end

    function jac_proj(v)
        u_curr = build_u(v)
        return ComplexF64[
            θ_shtiff(
                u_curr,
                p;
                char_shtiff=char_k,
                derivs=[[(i == global_index) ? 1 : 0 for i in 1:s]],
                faster=use_faster
            )
            for char_k in shtiff_chars, global_index in act_vars
        ]
    end

    if act_init === nothing
        v0 = ComplexF64.(collect(u_start[act_vars]))
    elseif (length(act_vars) == 1) && !(act_init isa AbstractVector)
        v0 = ComplexF64[act_init]
    else
        v0 = ComplexF64.(collect(act_init))
        if length(v0) != length(act_vars)
            error("project_to_image: длина act_init должна совпадать с длиной act_vars.")
        end
    end

    v_proj, diff_abs, iter, diff, diff_val = newton(f_proj, jac_proj, v0; tol=tol, maxiter=maxiter)

    #println("diff = ", diff)
    #println("diff_abs = ", diff_abs)

    #if diff_val > tol
    #    @warn "project_to_image: метод Ньютона мог не сойтись." diff_abs iter diff diff_val
    #end

    # !!! return build_u(v_proj)
    #??? Малый шаг Ньютона сам по себе не означает, что условие theta=0 выполнено.
    if !all(isfinite, v_proj) || !isfinite(diff_val) || diff_val > tol #???
        throw(DomainError(diff_val, "project_to_image: невязка theta остаётся больше tol")) #???
    end #???
    return build_u(v_proj) #???
end

function path_on_jacobian(u_11, u_12, u_22, p; tol=1e-10, step_num=100, frac_for_tang=100000)
    u2_dir = u_22 - u_12
    u2_step = u2_dir / step_num
    u_21 = u_11

    for i in 1:step_num

        u_21_prev = u_21
        norm_vect = [θ_shtiff([u_11, u_12], p; faster=use_faster, derivs=[[1, 0]], char_shtiff=[3, 5]),
            θ_shtiff([u_11, u_12], p; faster=use_faster, derivs=[[0, 1]], char_shtiff=[3, 5])]
        tangent_vect = [-norm_vect[2], norm_vect[1]]

        if abs(tangent_vect[1]/tangent_vect[2]) > frac_for_tang
            tangent_vect = [0, 1]
        end
        tangent_vect = tangent_vect / tangent_vect[2] * u2_step

        if norm(u2_step) < tol
            break
        end
        f = x -> θ_shtiff([x, u_12 + tangent_vect[2]], p; faster=use_faster, char_shtiff=[3, 5])
        dir_f = x -> θ_shtiff([x, u_12 + tangent_vect[2]], p; faster=use_faster, char_shtiff=[3, 5], derivs=[[1, 0]])

        time_in_newton = time()
        u_21, diff, iter = newton(f, dir_f, u_11 + tangent_vect[1]; tol=tol, maxiter=1000)
        time_in_newton = time() - time_in_newton

        if (diff > tol) || (abs((u_21-u_21_prev)) > 2*abs(norm(tangent_vect)))
            u_21 = path_on_jacobian(u_11, u_12, u_12 + u2_step, p, step_num=3, frac_for_tang=100000)
        end

        u_11 = u_21
        u_12 = u_12 + tangent_vect[2]
        #println("val = ", θ_shtiff(ComplexF64.([u_11, u_12]), p, char_shtiff = [3, 5]))
    end

    return u_21
end

function newton_method(u2, p, starting_point_index, step_num, newton_steps_count)
    u_curr = special_points(p)[starting_point_index]
    h2 = (u2 - u_curr[2]) / step_num   # шаг по второй координате на якобиане
    for i in 1:step_num
        u_curr[2] += h2
        for _ in 1:newton_steps_count
            u_curr[1] -= θ_shtiff(u_curr, p; char_shtiff=[3, 5], faster=use_faster) /
                         θ_shtiff(u_curr, p; derivs=[[1, 0]], char_shtiff=[3, 5], faster=use_faster)

        end
    end

    return u_curr
end

function decide_starting_point_particular(i)
    if i == M
        return 2
    end
    if i < M
        return 6
    end
    if i > M
        return 5
    end
end

const X_ON_SURFACE_LAST_PREFACTOR = Ref{Any}(nothing) #??? Только последний период: сетка почти не повторяет старые периоды.
const X_ON_SURFACE_PREFACTOR_LOCK = ReentrantLock() #???

function x_on_surface_prefactor(p, zero_index, one_index, infinity_index, added_index) #???
    key = p isa AbstractMatrix ? #???
        (size(p), Tuple(ComplexF64.(vec(p))), zero_index, one_index, infinity_index, added_index, USE_FASTER_THETA, THETA_EPS) : #???
        nothing #???
    if key !== nothing #???
        cached = lock(X_ON_SURFACE_PREFACTOR_LOCK) do #???
            entry = X_ON_SURFACE_LAST_PREFACTOR[] #???
            entry !== nothing && entry.key == key ? entry.value : nothing #???
        end #???
        cached !== nothing && return cached #???
    end #???
    th3_const = θ_shtiff([0, 0], p; char_shtiff=[infinity_index, added_index, 3, 5, one_index]) #???
    th4_const = θ_shtiff([0, 0], p; char_shtiff=[zero_index, added_index, 3, 5, one_index]) #???
    value = (th3_const / th4_const)^2 #???
    if key !== nothing #???
        lock(X_ON_SURFACE_PREFACTOR_LOCK) do #???
            X_ON_SURFACE_LAST_PREFACTOR[] = (; key, value) #???
        end #???
    end #???
    return value #???
end #???

include("surface_theta_log.jl") #??? Специальное вычисление для преобразованной матрицы периода.

function x_on_surface_new(u, p; zero_index=2, one_index=3, infinity_index=6, added_index=1, sec_added_index=1)
    τ = period_matrix(p) #???
    A, B, C, D = choose_modular_blocks(τ) #???
    if !is_identity_blocks(A, B, C, D) #???
        return _surface_x_with_cancelled_exponents(u, τ, A, B, C, D; #???
            zero_index=zero_index, one_index=one_index, infinity_index=infinity_index, #???
            added_index=added_index, sec_added_index=sec_added_index) #???
    end #???

    #str1 = zero_index * added_index * "35"
    #str2 = infinity_index * added_index * "35"
    #str3_const = str2 * one_index
    #str4_const = str1 * one_index

    th1 = θ_shtiff(u, p; char_shtiff=[zero_index, added_index, 3, 5, sec_added_index]) #, e_index=str1
    th2 = θ_shtiff(u, p; char_shtiff=[infinity_index, added_index, 3, 5, sec_added_index]) #, e_index=str2
    # !!! th3_const = θ_shtiff([0, 0], p; char_shtiff=[infinity_index, added_index, 3, 5, one_index])
    # !!! th4_const = θ_shtiff([0, 0], p; char_shtiff=[zero_index, added_index, 3, 5, one_index])
    prefactor = x_on_surface_prefactor(p, zero_index, one_index, infinity_index, added_index) #???

    # -1  [14352] [3435] k = 4, l = 1, j = 2, s = 3
    #     [34352] [1435]. ε(j) · (ε′(s) + ε′(l)),

    # !!! return -(th3_const / th4_const)^2 * (th1 / th2)^2
    return -prefactor * (th1 / th2)^2 #??? Отлично от формулы АБ; значение то же, константы не пересчитываются.
end

function fraction(zero_points, infinity_points, coefficient)
    # !!! P = coefficient * fromroots(zero_points)
    # !!! Q = fromroots(infinity_points)
    # !!! return z -> P(z) / Q(z)
    paired = min(length(zero_points), length(infinity_points)) #???
    return z -> begin #???
        value = coefficient
        for i in 1:paired
            value *= (z - zero_points[i]) / (z - infinity_points[i])
        end
        for i in (paired + 1):length(zero_points)
            value *= z - zero_points[i]
        end
        for i in (paired + 1):length(infinity_points)
            value /= z - infinity_points[i]
        end
        value
    end
end

#??? Дробно-линейное преобразование: a ↦ ∞, b ↦ -1, c ↦ 0.
function mobius_transform(a, b, c) #???
    all(isfinite, (a, b, c)) || throw(ArgumentError("a, b, c должны быть конечными числами")) #???
    (a != b && a != c && b != c) || throw(ArgumentError("a, b, c должны быть различными")) #???
    scale = (a - b) / (b - c) #???
    # !!! return x -> (x == a ? Inf : scale * (x - c) / (x - a))
    return x -> (x == a ? Inf : (isinf(x) ? scale : scale * (x - c) / (x - a))) #???
end #???

""" #???
    symmetric_mobius_normalization(x1, x4, x5, x6) #???

Вещественное дробно-линейное преобразование `T` с `T(x1)=-1`, #???
`T(x4)=1` и `T(x5)=-T(x6)`. Положительная ветвь масштаба не имеет #???
полюса между `x1` и `x4`. Возвращает `transform`, `inverse`, `gamma`. #???
""" #???
function symmetric_mobius_normalization(x1::Real, x4::Real, x5::Real, x6::Real) #???
    anchors = promote(float(x1), float(x4), float(x5), float(x6)) #???
    a, b, c, d = anchors #???
    all(isfinite, anchors) && length(unique(anchors)) == 4 || #???
        throw(ArgumentError("опорные точки x₁,x₄,x₅,x₆ должны быть конечными и различными: $anchors")) #???
    w5 = (c - a) / (b - c) #???
    w6 = (d - a) / (b - d) #???
    isfinite(w5) && isfinite(w6) && w5 != 0 && w6 != 0 && signbit(w5) == signbit(w6) || #???
        throw(ArgumentError("нет вещественного преобразования с T(x₅)=-T(x₆): точки=$anchors, w₅=$w5, w₆=$w6")) #???
    scale = sqrt(abs(w5)) * sqrt(abs(w6)) #??? Из симметрии следует scale²=w₅w₆.
    isfinite(scale) && scale > 0 || throw(ArgumentError("масштаб преобразования выродился: точки=$anchors")) #???
    transform(s) = ((s - a) - scale * (b - s)) / ((s - a) + scale * (b - s)) #???
    inverse(u) = ((1 - u) * a + scale * (1 + u) * b) / ((1 - u) + scale * (1 + u)) #???
    image5, image6 = transform(c), transform(d) #???
    all(isfinite, (image5, image6)) && image5 != 0 || #???
        throw(ArgumentError("образы x₅,x₆ не конечны либо 1/γ вырождается: точки=$anchors")) #???
    abs(image5 + image6) <= 1e-10 * max(1, abs(image5), abs(image6)) || #???
        throw(ArgumentError("не удалось получить симметричные образы x₅,x₆: $((image5, image6))")) #???
    # !!! return (; transform, inverse, gamma=inv(image5), scale)
    pole = scale == 1 ? Inf : (a - scale * b) / (1 - scale) #???
    return (; transform, inverse, gamma=inv(image5), scale, pole) #???
end #???

""" #???
    six_point_mobius_normalization(x1, x2, x3, x4, x5, x6) #???

Нормировка с `T(x3)=1`, `T(x6)=-1`, `T(x1)=-T(x2)`. #???
Из двух вещественных ветвей выбирается та, для которой #???
`α=T(x5)` и `β=T(x4)` лежат строго между -1 и 1. #???
""" #???
function six_point_mobius_normalization(x1::Real, x2::Real, x3::Real, #???
                                         x4::Real, x5::Real, x6::Real) #???
    anchors = promote(float(x1), float(x2), float(x3), #???
                      float(x4), float(x5), float(x6)) #???
    all(isfinite, anchors) && length(unique(anchors)) == 6 || #???
        throw(ArgumentError("для нормировки нужны шесть конечных различных точек: $anchors")) #???
    a, b = anchors[3], anchors[6] #??? Эти точки переходят в 1 и -1.
    base = symmetric_mobius_normalization(a, b, anchors[1], anchors[2]) #??? Модуль масштаба задаёт симметрия x₁,x₂.
    candidates = NamedTuple[] #???
    for scale in (base.scale, -base.scale) #??? Проверяем обе вещественные ветви.
        transform(s) = -((s - a) - scale * (b - s)) / #???
                       ((s - a) + scale * (b - s)) #???
        inverse(u) = ((1 + u) * a + scale * (1 - u) * b) / #???
                     ((1 + u) + scale * (1 - u)) #???
        images = transform.(anchors) #???
        all(isfinite, images) || continue #???
        abs(images[1] + images[2]) <= 1e-10 * max(1, abs(images[1]), abs(images[2])) || continue #???
        images[1] != 0 || continue #???
        α, β = images[5], images[4] #???
        pole = scale == 1 ? Inf : (a - scale * b) / (1 - scale) #???
        push!(candidates, (; transform, inverse, alpha=α, beta=β, #???
                            gamma=inv(images[1]), scale, pole, images)) #???
    end #???
    valid = filter(c -> abs(c.alpha) < 1 && abs(c.beta) < 1, candidates) #???
    isempty(valid) && throw(ArgumentError( #???
        "нет вещественной ветви с |α|<1 и |β|<1 для точек $anchors; образы кандидатов: $(getproperty.(candidates, :images))")) #???
    return first(valid) #???
end #???

function fraction_der(zero_points, infinity_points, coefficient)
    P = coefficient * fromroots(zero_points)
    der_P = derivative(P)
    Q = fromroots(infinity_points)
    der_Q = derivative(Q)
    return z -> (der_P(z) * Q(z) - P(z) * der_Q(z)) / Q(z)^2
end


"""
u_ij, i -- номер точки, j -- номер координаты в точке
"""
function path_on_jacobian_new(u_11, u_12, u_21, p; tol=1e-10, step_num=100, frac_for_tang=100000)
    u1_dir = u_21 - u_11
    u1_step = u1_dir / step_num
    u_21_curr = u_11
    u_22_curr = u_12
    for i in 1:step_num
        u_21_prev = u_21_curr
        u_22_prev = u_22_curr

        norm_vect = [θ_shtiff([u_21_curr, u_22_curr], p; faster=use_faster, derivs=[[1, 0]], char_shtiff=[3, 5]),
            θ_shtiff([u_21_curr, u_22_curr], p; faster=use_faster, derivs=[[0, 1]], char_shtiff=[3, 5])]
        tangent_vect = [-norm_vect[2], norm_vect[1]]

        if abs(tangent_vect[2]/tangent_vect[1]) > frac_for_tang
            tangent_vect = [1, 0]
        end
        tangent_vect = tangent_vect / tangent_vect[1] * u1_step

        if norm(u1_step) < tol
            break
        end

        u_trial = [u_21_curr, u_22_curr] .+ tangent_vect
        u_proj = project_to_image(u_trial, p; shtiff_chars=[[3, 5]], act_vars=[2], tol=tol, maxiter=100)

        residual = norm(θ_shtiff(u_proj, p; char_shtiff=[3, 5], faster=use_faster))
        if (residual > tol) #|| (abs(u_proj[2] - u_22_prev) > 2 * abs(norm(tangent_vect)))
            u_22_curr = path_on_jacobian_new(u_21_prev, u_22_prev, u_21_prev + u1_step, p; tol=tol, step_num=3, frac_for_tang=frac_for_tang)
            u_21_curr = u_21_prev + u1_step
            continue
        end

        u_21_curr = u_proj[1]
        u_22_curr = u_proj[2]
    end

    return u_22_curr
end

function shtiff_fraction_new(t, n, p, m; tol=1e-10)
    points = ComplexF64[k * t * im for k in 0:n] #* (-0.25) + 0.25
    points_count = n + 1

    for i in 1:points_count
        points[i] -= 1
        points[i] /= -4
    end


    zero_points1 = []
    infinity_points1 = []

    #u_1_curr = path_on_jacobian(0, 0, 0.5, RiemannMatrix(p, siegel=false))
    u_2_curr = 0
    u_curr = [0., 0.]

    # !!! u_2_curr_down = path_on_jacobian_new(u_curr[1], u_curr[2], 0.25*im*n*t, p)
    u_2_curr_down = path_on_jacobian_new(u_curr[1], u_curr[2], -0.25*im*n*t, p) #???
    u_curr_down = [- 0.25*im*n*t, u_2_curr_down]

    for i in 1:points_count
        # начальная точка, от которой идем методом ньютона до нашей
        #starting_point_index = decide_starting_point_particular(i - 1)

        # решение уравнения theta[35](u) = 0 с u2 из points
        time_in_newton = time()
        #u = newton_method(points[i], p, starting_point_index, newton_number, 5)
        u_22 = path_on_jacobian_new(u_curr[1], u_curr[2], points[i], p)
        time_in_newton = time() - time_in_newton
        println("time in newton ", time_in_newton)

        # образ при отображении КШ(u2)
        x = x_on_surface_new([points[i], u_22], p)
        #u_curr = [points[i], u_22]
        println("u_curr = ", [points[i], u_22])
        println("θ_shtiff(u_curr, p; char_shtiff = [3, 5]) = ", θ_shtiff([points[i], u_22], p; char_shtiff = [3, 5]))

        #println("here")

        # образы точек i=0,n лежат на вещественной прямой, так что это уже нули или полюса, сопряженных к ним нет
        if i == 1
            push!(zero_points1, x)

        elseif i == n + 1
            if (i - 1) % 2 == 0
                push!(zero_points1, x)
            else
                push!(infinity_points1, x)
            end

            # if разрез длиннее 0.25, то обе точки x, x_added должны лежать на вещественной прямой, так что их не сопрягаем
            elseif i == m + 1
            # поскольку(в будущей версии кода)
            # находимся в случае, когда разрез длинный, i=m соответствует два вещественных образа при КШ
                u_22_down = path_on_jacobian_new(u_curr_down[1], u_curr_down[2], points[i], p)
                u_down = [points[i], u_22_down]
                #u_curr_down = copy(u_down)
                x_added = x_on_surface_new(u_down, p)
                if abs(imag(x)) > tol && abs(x_added - x) < abs(x_added - conj(x)) #???
                    x_added = conj(x) #??? На нижнем пути найден тот же комплексный образ вместо сопряжённого.
                end

            # для спокойствия посмотрим, что точки вещественны
            # println(x, " ", x_added)

            # добавляем оба вещественных образа в ноль или полюс понятным образом
            #println()
            #println("abs(x_added - x) = ", abs(x_added - x))
            #println("x_added = ", x_added, " x = ", x)
            #println()
            if (i - 1) % 2 == 0
                # println("two zeros")
                if abs(x_added - x) < tol
                    push!(zero_points1, x)
                    push!(zero_points1, conj(x_added))
                else
                    push!(zero_points1, x)
                    push!(zero_points1, x_added)
                end
            else
                if abs(x_added - x) < tol
                    push!(infinity_points1, x)
                    push!(infinity_points1, conj(x_added))
                else
                    push!(infinity_points1, x)
                    push!(infinity_points1, x_added)
                end
                #push!(infinity_points1, x)
                #push!(infinity_points1, x_added)
            end

            # if разрез меньше 0.25, то там надо добавить сопряженные тоже
            # код здесь будет такой, как ниже для остальных точек i != 0, m, n ......
            # .....
        elseif i < m + 1
            # поскольку тут осталось учесть точки, заведомо имеющие комплексный образ, надо добавить сопряженные к ним.
            if (i - 1) % 2 == 0
                push!(zero_points1, x)
                push!(zero_points1, conj(x))
            else
                push!(infinity_points1, x)
                push!(infinity_points1, conj(x))
            end
        else
            u_22_down = path_on_jacobian_new(u_curr_down[1], u_curr_down[2], points[i], p)
            u_down = [points[i], u_22_down]
            #u_curr_down = copy(u_down)
            x_added = x_on_surface_new(u_down, p)

            if (i - 1) % 2 == 0
                push!(zero_points1, x_added)
                push!(zero_points1, conj(x_added))
            else
                push!(infinity_points1, x_added)
                push!(infinity_points1, conj(x_added))
            end

        end


        # println(i)
    end

    # здесь можно преобразовать нули и полюса с помощью дробно-линейного преобразования
    zero_points = ComplexF64.(zero_points1)
    infinity_points = ComplexF64.(infinity_points1)
    frac = fraction(zero_points, infinity_points, 1)

    # нормировка значения дроби: дробь(1) = 1
    coefficient = 1 / frac(1)[1]

    # "существенная" область определения дроби -- отбросим самые крайние нули и полюса
    x_max_z = sort(abs.(zero_points))
    x_max_p = sort(abs.(infinity_points))

    # полная информация о дроби
    # println("zero_points, infinity_points = $(zero_points, infinity_points), coefficient = $coefficient")

    return zero_points, infinity_points, coefficient, max(x_max_z[end], x_max_p[end])
end

"""
"""

function shtiff_fraction(t, n, p, m; tol=1e-10)
    points = ComplexF64[k * t * im for k in 0:n] #* (-0.25) + 0.25
    points_count = n + 1

    for i in 1:points_count
        points[i] -= 1
        points[i] /= -4
    end


    zero_points1 = []
    infinity_points1 = []

    u_1_curr = path_on_jacobian(0, 0, 0.5, p)
    u_curr = [u_1_curr, 0.5]

    u_1_curr_down = path_on_jacobian(u_curr[1], u_curr[2], 0.5 - 0.25*im*n*t, p)
    u_curr_down = [u_1_curr_down, 0.5 - 0.25*im*n*t]

    for i in 1:points_count
        # начальная точка, от которой идем методом ньютона до нашей
        #starting_point_index = decide_starting_point_particular(i - 1)

        # решение уравнения theta[35](u) = 0 с u2 из points
        time_in_newton = time()
        #u = newton_method(points[i], p, starting_point_index, newton_number, 5)
        u_21 = path_on_jacobian(u_curr[1], u_curr[2], points[i], p)
        time_in_newton = time() - time_in_newton
        println("time in newton ", time_in_newton)

        # образ при отображении КШ(u2)
        x = x_on_surface_new([u_21, points[i]], p)
        u_curr = [u_21, points[i]]

        println("u_curr = ", u_curr)
        println("θ_shtiff(u_curr, p; char_shtiff = [3, 5]) = ", θ_shtiff(u_curr, p; char_shtiff = [3, 5]))
        #println("here")

        # образы точек i=0,n лежат на вещественной прямой, так что это уже нули или полюса, сопряженных к ним нет
        if i == 1
            push!(zero_points1, x)

        elseif i == n + 1
            if (i - 1) % 2 == 0
                push!(zero_points1, x)
            else
                push!(infinity_points1, x)
            end

            # if разрез длиннее 0.25, то обе точки x, x_added должны лежать на вещественной прямой, так что их не сопрягаем
        elseif i == m + 1
            # поскольку(в будущей версии кода)
            # находимся в случае, когда разрез длинный, i=m соответствует два вещественных образа при КШ
            u_21_down = path_on_jacobian(u_curr_down[1], u_curr_down[2], points[i], p)
            u_down = [u_21_down, points[i]]
            u_curr_down = copy(u_down)
            x_added = x_on_surface_new(u_down, p)

            # для спокойствия посмотрим, что точки вещественны
            # println(x, " ", x_added)

            # добавляем оба вещественных образа в ноль или полюс понятным образом
            #println()
            #println("abs(x_added - x) = ", abs(x_added - x))
            #println("x_added = ", x_added, " x = ", x)
            #println()
            if (i - 1) % 2 == 0
                # println("two zeros")
                if abs(x_added - x) < tol
                    push!(zero_points1, x)
                    push!(zero_points1, conj(x_added))
                else
                    push!(zero_points1, x)
                    push!(zero_points1, x_added)
                end
            else
                if abs(x_added - x) < tol
                    push!(infinity_points1, x)
                    push!(infinity_points1, conj(x_added))
                else
                    push!(infinity_points1, x)
                    push!(infinity_points1, x_added)
                end
                #push!(infinity_points1, x)
                #push!(infinity_points1, x_added)
            end

            # if разрез меньше 0.25, то там надо добавить сопряженные тоже
            # код здесь будет такой, как ниже для остальных точек i != 0, m, n ......
            # .....
        elseif i < m + 1
            # поскольку тут осталось учесть точки, заведомо имеющие комплексный образ, надо добавить сопряженные к ним.
            if (i - 1) % 2 == 0
                push!(zero_points1, x)
                push!(zero_points1, conj(x))
            else
                push!(infinity_points1, x)
                push!(infinity_points1, conj(x))
            end
        else
            u_21_down = path_on_jacobian(u_curr_down[1], u_curr_down[2], points[i], p)
            u_down = [u_21_down, points[i]]
            u_curr_down = copy(u_down)
            x_added = x_on_surface_new(u_down, p)

            if (i - 1) % 2 == 0
                push!(zero_points1, x_added)
                push!(zero_points1, conj(x_added))
            else
                push!(infinity_points1, x_added)
                push!(infinity_points1, conj(x_added))
            end

        end


        # println(i)
    end

    # здесь можно преобразовать нули и полюса с помощью дробно-линейного преобразования
    zero_points = ComplexF64.(zero_points1)
    infinity_points = ComplexF64.(infinity_points1)
    frac = fraction(zero_points, infinity_points, 1)

    # нормировка значения дроби: дробь(1) = 1
    coefficient = 1 / frac(1)[1]

    # "существенная" область определения дроби -- отбросим самые крайние нули и полюса
    x_max_z = sort(abs.(zero_points))
    x_max_p = sort(abs.(infinity_points))

    # полная информация о дроби
    # println("zero_points, infinity_points = $(zero_points, infinity_points), coefficient = $coefficient")

    return zero_points, infinity_points, coefficient, max(x_max_z[end], x_max_p[end])
end

function shtiff_fraction_gen_1(t, n, p; tol=1e-10)
    n isa Integer && n >= 1 || throw(ArgumentError("n должно быть положительным целым числом")) #???
    # !!! coeff = -im * θ([0], p)/θ([im * t/2], p)
    # !!! func = x -> coeff * exp(pi * im * x/2) * θ([(x + im * t +1)/2], p)/θ([(x + 1)/2], p)
    τ_total = period_matrix(p)[1, 1] #??? Период соответствует высоте всего прямоугольника.
    isapprox(τ_total, im * n * t; atol=tol, rtol=tol) || throw(ArgumentError("Для точек k*im*t период p[1,1] должен быть im*n*t")) #???
    coeff = -im * θ([0], p) / θ([τ_total / 2], p) #???
    func = x -> coeff * exp(pi * im * x / 2) * θ([(x + τ_total + 1) / 2], p) / θ([(x + 1) / 2], p) #???
    x1_value = func(-1 + τ_total) #???
    x2_value = func(-1) #???
    x3_value = func(1) #???
    for boundary_value in (x1_value, x2_value, x3_value) #???
        isfinite(boundary_value) || throw(DomainError(boundary_value, "Образ граничной точки должен быть конечным")) #???
        if abs(imag(boundary_value)) > tol * max(1.0, abs(real(boundary_value))) #???
            throw(DomainError(boundary_value, "Образ граничной точки должен быть действительным")) #???
        end #???
    end #???
    x1, x2, x3 = real(x1_value), real(x2_value), real(x3_value) #???
    T = mobius_transform(x1, x2, x3) #??? Граничные точки переходят в ∞, -1 и 0.
    #= !!! Старое формирование списков: push! добавлял вектор как один элемент,
            а при p[1,1] = im*n*t верхняя точка im*n*t является полюсом func.
    if n % 2 == 0
        zeroes = func.(ComplexF64[k * t * 2 * im for k in 0:n/2]) #* (-0.25) + 0.25
        infinite = func.(ComplexF64[k * t * 2 * im + t * im for k in 0:n/2-1]) #* (-0.25) + 0.25
        
        push!(zeroes, conj.(zeroes))
        push!(infinite, conj.(infinite))
        return zeroes, infinite
    else
        zeroes = func.(ComplexF64[k * t * 2 * im for k in -n/2+1:n/2]) #* (-0.25) + 0.25
        infinite = func.(ComplexF64[k * t * 2 * im + t * im for k in -n/2:n/2-1]) #* (-0.25) + 0.25
        return zeroes, infinite
    end
    =#
    #= !!! До замены координаты через T верхняя вершина не была конечным корнем.
    #??? При p = im*n*t верхняя вершина переходит в ∞: её нуль или полюс задаёт разность степеней.
    zeroes = ComplexF64[func(im * t * k) for k in 0:2:(n - 1)] #???
    infinite = ComplexF64[func(im * t * k) for k in 1:2:(n - 1)] #???
    # !!! append!(zeroes, conj.(zeroes))
    zeroes[1] = 0.0 + 0.0im #??? Граничный образ func(0) = 0 учитывается один раз.
    append!(zeroes, conj.(zeroes[2:end])) #??? Сопряжённые добавляем только для внутренних точек.
    append!(infinite, conj.(infinite)) #???
    return zeroes, infinite #???
    =#
    zeroes = ComplexF64[] #???
    infinite = ComplexF64[] #???
    for k in 0:n #???
        #??? На концах используем точные значения func(0)=0 и func(τ_total)=∞.
        mapped = k == 0 ? T(0.0) : (k == n ? T(Inf) : T(func(im * t * k))) #???
        isfinite(mapped) || throw(DomainError(mapped, "Образ точки k=$k не конечен")) #???
        target = iseven(k) ? zeroes : infinite #???
        push!(target, mapped) #???
        if k != 0 && k != n #???
            push!(target, conj(mapped)) #???
        end #???
    end #???
    return zeroes, infinite #???
end

function draw_shtiff(zero_points, infinity_points, coefficient, x_max)
    @eval using Plots
    x = range(-x_max, x_max; length=100_000)
    frac = fraction(zero_points, infinity_points, coefficient)
    y = real.(frac.(x))

    p = plot(
        x, y,
        xlabel="x",
        ylabel="Re(f(x))",
        label="",
        legend=false
    )

    display(p)
    return p
end

#=
function eversion_rectangular(x, τ, p)
    return -im * exp(-π * im * τ / 4) *
           θ(0.0 + 0.0im, p) / θ(τ/2, p) *
           θ(x/2 + 1/2 + τ/2, p) / θ(x/2 + τ/2, p)
end
=#

"""
Строит график y = f(x) в отдельном интерактивном окне.

Масштаб:
- колесо мыши — приблизить/отдалить
- левая кнопка + протянуть — зум в прямоугольник
- правая кнопка + тянуть — перемещение
- Ctrl + левая кнопка — сброс масштаба
"""
    
# !!! using GLMakie
#??? Для сохранения 2D-графика без оконной системы можно выбрать CairoMakie.
if lowercase(get(ENV, "SHTIFFEL_PLOT_BACKEND", "gl")) == "cairo" #???
    using CairoMakie #???
    CairoMakie.activate!() #???
else #???
    using GLMakie #???
end #???
import Makie

#= !!! Прежний график только функции f, без кусочно-постоянной функции.
function draw_graph(f, x_min, x_max; n=1000)
    x_values = range(x_min, x_max; length=n)
    y_values = f.(x_values)

    fig = Makie.Figure()
    ax = Makie.Axis(fig[1, 1])

    Makie.lines!(ax, x_values, y_values)
    Makie.xlims!(ax, x_min, x_max) #???

    display(fig)
    return fig
end
=#

# !!! intervals и values задают три участка кусочно-постоянной функции.
#??? intervals и values задают три или четыре участка кусочно-постоянной функции.
#??? Бесконечные концы интервалов обрезаются видимыми пределами x_min, x_max.
#??? Концы каждого интервала можно передавать в любом порядке.
#??? Разрываем линию между соседними точками по разные стороны полюса.
function _graph_line_with_gaps(x_values, y_values) #???
    ys = map(y_values) do value #???
        value isa Real || throw(ArgumentError("f(x) должна возвращать действительные значения")) #???
        number = Float64(value) #???
        isfinite(number) ? number : NaN #???
    end #???
    sample_step = max(1, cld(length(ys), 2048)) #???
    sample = sort!([abs(ys[i]) for i in 1:sample_step:length(ys) if isfinite(ys[i])]) #???
    typical_height = isempty(sample) ? 0.0 : sample[cld(length(sample), 2)] #???
    gaps = falses(length(ys) - 1) #???
    for i in eachindex(gaps) #???
        left, right = ys[i], ys[i + 1] #???
        isfinite(left) && isfinite(right) && signbit(left) != signbit(right) || continue #???
        jump = abs(right - left) #???
        before = i > 1 && isfinite(ys[i - 1]) ? abs(left - ys[i - 1]) : 0.0 #???
        after = i + 2 <= length(ys) && isfinite(ys[i + 2]) ? abs(ys[i + 2] - right) : 0.0 #???
        gaps[i] = jump > 6 * typical_height || #???
                  (before > 0 && after > 0 && jump > 2.5 * max(before, after)) #???
    end #???
    any(gaps) || return x_values, ys #???
    xs_plot = Float64[] #???
    ys_plot = Float64[] #???
    sizehint!(xs_plot, length(ys) + count(gaps)) #???
    sizehint!(ys_plot, length(ys) + count(gaps)) #???
    for i in eachindex(ys) #???
        push!(xs_plot, Float64(x_values[i])) #???
        push!(ys_plot, ys[i]) #???
        if i < length(ys) && gaps[i] #???
            push!(xs_plot, (Float64(x_values[i]) + Float64(x_values[i + 1])) / 2) #???
            push!(ys_plot, NaN) #???
        end #???
    end #???
    return xs_plot, ys_plot #???
end #???

# !!! function draw_graph(f, x_min, x_max; n=1000, intervals=nothing, values=nothing, ylims=nothing)
function draw_graph(f, x_min, x_max; n=1000, intervals=nothing, values=nothing, #???
                    ylims=nothing, wrapped_intervals=()) #???
    isfinite(x_min) && isfinite(x_max) && x_min < x_max || #???
        throw(ArgumentError("нужны конечные x_min < x_max")) #???
    n isa Integer && n >= 2 || throw(ArgumentError("n должно быть целым числом не меньше 2")) #???
    if ylims !== nothing #???
        ylims isa Tuple && length(ylims) == 2 && #???
            all(y -> y isa Real && isfinite(y), ylims) && ylims[1] < ylims[2] || #???
            throw(ArgumentError("ylims должно иметь вид (y_min, y_max), где y_min < y_max — конечные числа")) #???
    end #???
    (intervals === nothing) == (values === nothing) || #???
        throw(ArgumentError("intervals и values нужно передать вместе")) #???
    #??? Участок, содержащий полюс преобразования, состоит из двух лучей на вещественной оси.
    wrapped_intervals isa Tuple && all(k -> k isa Integer, wrapped_intervals) || #???
        throw(ArgumentError("wrapped_intervals должно быть кортежем номеров интервалов")) #???
    intervals === nothing && !isempty(wrapped_intervals) && #???
        throw(ArgumentError("wrapped_intervals требует intervals")) #???
    if intervals !== nothing #???
        # !!! length(intervals) == 3 && length(values) == 3 ||
        # !!!     throw(ArgumentError("нужно ровно три интервала и три значения"))
        length(intervals) in (3, 4) && length(values) == length(intervals) || #???
            throw(ArgumentError("нужно три или четыре интервала и столько же значений")) #???
        all(k -> 1 <= k <= length(intervals), wrapped_intervals) || #???
            throw(ArgumentError("номер в wrapped_intervals выходит за число интервалов")) #???
        for (k, interval) in enumerate(intervals) #???
            length(interval) == 2 || throw(ArgumentError("интервал $k должен иметь два конца")) #???
            a, b = interval #???
            # !!! a isa Real && b isa Real && !isnan(a) && !isnan(b) && a < b ||
            # !!!     throw(ArgumentError("интервал $k должен иметь вид (a, b), где a < b"))
            # !!! a isa Real && b isa Real && !isnan(a) && !isnan(b) && a != b ||
            # !!!     throw(ArgumentError("интервал $k должен иметь два разных действительных конца"))
            a isa Real && b isa Real && !isnan(a) && !isnan(b) || #???
                throw(ArgumentError("интервал $k должен иметь два действительных конца, не NaN")) #???
            a == b && @warn "Интервал нулевой длины пропущен" номер=k граница=a #???
            values[k] isa Real && isfinite(values[k]) || #???
                throw(ArgumentError("значение $k должно быть конечным действительным числом")) #???
        end #???
    end #???

    x_values = range(x_min, x_max; length=n) #???
    y_values = f.(x_values) #???
    x_line, y_line = _graph_line_with_gaps(x_values, y_values) #???
    fig = Makie.Figure() #???
    ax = Makie.Axis(fig[1, 1]) #???
    # !!! Makie.lines!(ax, x_values, y_values)
    Makie.lines!(ax, x_line, y_line) #???
    if intervals !== nothing #???
        for (k, interval) in enumerate(intervals) #???
            a, b = interval #???
            # !!! left, right = max(x_min, min(a, b)), min(x_max, max(a, b))
            pieces = k in wrapped_intervals ? #???
                     ((x_min, min(a, b)), (max(a, b), x_max)) : #???
                     ((min(a, b), max(a, b)),) #???
            # !!! color = k == 1 ? :red : :green
            color = (k == 1 || (length(intervals) == 4 && k == 4)) ? :red : :green #???
            for (piece_left, piece_right) in pieces #???
                left, right = max(x_min, piece_left), min(x_max, piece_right) #???
                left < right || continue #???
                Makie.lines!(ax, [left, right], [values[k], values[k]]; #???
                             color=color, linewidth=3) #???
            end #???
        end #???
    end #???
    Makie.xlims!(ax, x_min, x_max) #???
    ylims === nothing || Makie.ylims!(ax, ylims...) #???
    display(fig) #???
    return fig #???
end #???

#??? Численный максимум действительной функции на конечном отрезке.
#??? Сетка находит несколько пиков, затем каждый уточняется золотым сечением.
"""
    find_maximum(f, x_min, x_max; n=10_001, iterations=40)

Возвращает `(x=..., value=...)` — приближённый максимум `f` на `[x_min, x_max]`.
`n` задаёт число точек первоначальной сетки. Для функции с полюсом или
разрывом на отрезке конечного максимума может не быть.
"""
function find_maximum(f, x_min::Real, x_max::Real; n::Integer=10_001,
                      iterations::Integer=40) #???
    isfinite(x_min) && isfinite(x_max) && x_min < x_max || #???
        throw(ArgumentError("нужны конечные x_min < x_max")) #???
    n >= 3 || throw(ArgumentError("n должно быть не меньше 3")) #???
    iterations >= 0 || throw(ArgumentError("iterations не может быть отрицательным")) #???

    function finite_value(x) #???
        y = f(x) #???
        y isa Real && isfinite(y) || #???
            throw(DomainError(x, "f(x) должна быть конечным действительным числом")) #???
        value = Float64(y) #???
        isfinite(value) || throw(DomainError(x, "f(x) выходит за диапазон Float64")) #???
        return value #???
    end #???

    xs = range(Float64(x_min), Float64(x_max); length=n) #???
    ys = [finite_value(x) for x in xs] #???
    best_index = argmax(ys) #???
    best_x, best_value = xs[best_index], ys[best_index] #???
    ratio = (sqrt(5.0) - 1.0) / 2.0 #???

    for i in 2:(n - 1) #???
        ys[i] >= ys[i - 1] && ys[i] >= ys[i + 1] && #???
            (ys[i] > ys[i - 1] || ys[i] > ys[i + 1]) || continue #???
        left, right = xs[i - 1], xs[i + 1] #???
        c, d = right - ratio * (right - left), left + ratio * (right - left) #???
        fc, fd = finite_value(c), finite_value(d) #???
        for _ in 1:iterations #???
            # !!! right - left <= 16 * eps(max(abs(left), abs(right), 1.0)) && break
            right - left <= 16 * Base.eps(max(abs(left), abs(right), 1.0)) && break #???
            if fc < fd #???
                left, c, fc = c, d, fd #???
                d = left + ratio * (right - left) #???
                fd = finite_value(d) #???
            else #???
                right, d, fd = d, c, fc #???
                c = right - ratio * (right - left) #???
                fc = finite_value(c) #???
            end #???
        end #???
        if fc > best_value #???
            best_x, best_value = c, fc #???
        end #???
        if fd > best_value #???
            best_x, best_value = d, fd #???
        end #???
    end #???
    return (x=best_x, value=best_value) #???
end #???

#print("x_on_surface_new ", x_on_surface_new(u, p))

#special_points(im * [3.141592653589793 2.71828182845904; 2.718281828459045 3.141592653589793])  # пример вызова функции для π
