using Theta
using LinearAlgebra

USE_FASTER_THETA = false
# !!! const THETA_EPS = 1.0e-15 #??? Общий допуск для Theta.jl.
const THETA_EPS = 5.0e-16 #??? Допуск усечения Theta.jl; точность Float64 остаётся ограничением.
# !!! const THETA_DERIV_STEP = 1.0e-3 #??? Шаг для пятиузловой формулы производной.
# !!! const THETA_DERIV_STEP = 1.0e-4 #??? Уменьшенный шаг для пятиузловой формулы производной.
const THETA_DERIV_STEP = 1.0e-6 #??? Шаг численных производных тета-функции.

"""
Возвращает два столбца матрицы тэта-характеристики.
"""
function e_columns(t)
    return t[:, 1], t[:, 2]
end

"""
Возвращает набор специальных тэта-характеристик, используемых в коде.
"""
function special_chars()
    e1 = [0. 0.; 0. 0.]
    e2 = [1. 0.; 0. 0.]
    e3 = [1. 1.; 0. 0.]
    e4 = [0. 1.; 1. 0.]
    e5 = [0. 1.; 1. 1.]
    e6 = [0. 1.; 0. 1.]
    return [e1, e2, e3, e4, e5, e6]
end

"""
Собирает матрицу характеристики из обычной матрицы char и списка специальных характеристик char_shtiff.
"""
function shtiff_char_matrix(g; char=nothing, char_shtiff=[])
    char_matrix = if char === nothing || length(char) == 0
        zeros(Int, g, g)
    else
        Int.(copy(char))
    end

    if length(char_shtiff) > 0
        spe_char = special_chars()
        for idx in char_shtiff
            char_matrix .+= Int.(spe_char[idx])
        end
    end

    return char_matrix
end

include("siegel.jl")

"""
Численно дифференцирует функцию по переменной u в направлении deriv.
"""
#= !!! Прежняя двухузловая формула с h=1e-6 усиливала ошибки округления значений theta.
function finite_diff_u(f, u, p, deriv; h=1e-6)
    direction = ComplexF64.(deriv)
    return (f(u .+ h .* direction, p) - f(u .- h .* direction, p)) / (2h)
end
=#
function finite_diff_u(f, u, p, deriv; h=THETA_DERIV_STEP) #???
    direction = ComplexF64.(deriv) #???
    step = h .* direction #???
    inner = f(u .+ step, p) - f(u .- step, p) #???
    outer = f(u .+ 2 .* step, p) - f(u .- 2 .* step, p) #???
    return (8 * inner - outer) / (12h) #???
end #???

"""
Численно дифференцирует функцию по матрице периода p в направлении элемента deriv_t = [i, j].
"""
#= !!! Прежняя двухузловая формула с h=1e-6 усиливала ошибки округления значений theta.
function finite_diff_p(f, u, p, deriv_t; h=1e-6)
    p_plus = copy(p)
    p_minus = copy(p)
    i, j = deriv_t
    p_plus[i, j] += h
    p_minus[i, j] -= h
    if i != j
        p_plus[j, i] += h
        p_minus[j, i] -= h
    end
    return (f(u, p_plus) - f(u, p_minus)) / (2h)
end
=#
function finite_diff_p(f, u, p, deriv_t; h=THETA_DERIV_STEP) #???
    i, j = deriv_t #???
    function shifted_value(shift) #???
        p_shifted = copy(p) #???
        p_shifted[i, j] += shift #???
        if i != j #???
            p_shifted[j, i] += shift #???
        end #???
        return f(u, p_shifted) #???
    end #???
    inner = shifted_value(h) - shifted_value(-h) #???
    outer = shifted_value(2h) - shifted_value(-2h) #???
    return (8 * inner - outer) / (12h) #???
end #???

"""
Рекурсивно применяет численные производные по u и p к функции f.
"""
# !!! function apply_transformed_derivs(f, u, p, derivs, derivs_t; h=1e-6)
function apply_transformed_derivs(f, u, p, derivs, derivs_t; h=THETA_DERIV_STEP) #???
    if length(derivs) > 0
        deriv = derivs[1]
        return finite_diff_u(
            (u_curr, p_curr) -> apply_transformed_derivs(f, u_curr, p_curr, derivs[2:end], derivs_t; h=h),
            u,
            p,
            deriv;
            h=h,
        )
    end
    if length(derivs_t) > 0
        deriv_t = derivs_t[1]
        return finite_diff_p(
            (u_curr, p_curr) -> apply_transformed_derivs(f, u_curr, p_curr, derivs, derivs_t[2:end]; h=h),
            u,
            p,
            deriv_t;
            h=h,
        )
    end
    return f(u, p)
end

"""
Вычисляется значение тета-функции в точке u при параметре p.

`derivs` задаёт производные по переменным u.
`derivs_t` задаёт производные по элементам матрицы периода.
`char` задаёт тэта-характеристику в виде двух столбцов.

Если `faster=false`, то сначала применяется проектное преобразование Зигеля,
а уже затем вызывается `RiemannMatrix(..., siegel=false, ϵ=...)`.
Это правило одинаково работает и для обычной матрицы `p`, и для `RiemannMatrix`.

Если `faster=true`, то используется прямой вызов библиотеки Theta.jl.
Этот режим действительно быстрый только когда `p` уже имеет тип `RiemannMatrix`.
"""
# !!! function θ(u, p; derivs=[], derivs_t=[], char=[], faster=USE_FASTER_THETA, ϵ=1.0e-12)
function θ(u, p; derivs=[], derivs_t=[], char=[], faster=USE_FASTER_THETA, ϵ=THETA_EPS) #???
    # !!! if faster
    if faster && length(u) > 1 #??? Для рода 1 используем обход ошибки в Theta.jl через матрицу 2×2.
        # !!! if p isa RiemannMatrix
        if p isa Theta.RiemannMatrix #???
            # !!! return theta(u, p, derivs=derivs, derivs_t=derivs_t, char=char)
            return Theta.theta(u, p, derivs=derivs, derivs_t=derivs_t, char=char) #???
        end
        τ = period_matrix(p)
        # !!! return theta(u, RiemannMatrix(τ, siegel=false, ϵ=ϵ), derivs=derivs, derivs_t=derivs_t, char=char)
        # !!! return Theta.theta(u, theta_matrix_cached(τ; ϵ=ϵ), derivs=derivs, derivs_t=derivs_t, char=char)
        derivative_order = length(derivs) + 2 * length(derivs_t) #??? По уравнению теплопроводности каждая производная по τ даёт две по u.
        return Theta.theta(u, theta_matrix_cached(τ; ϵ=ϵ, nderivs=derivative_order), #???
                           derivs=derivs, derivs_t=derivs_t, char=char) #???
    end

    if isempty(derivs) && isempty(derivs_t)
        return theta_value_auto(u, p; char=char, ϵ=ϵ)
    end

    τ = period_matrix(p)
    value_function = (u_curr, p_curr) -> theta_value_auto(u_curr, p_curr; char=char, ϵ=ϵ)
    return apply_transformed_derivs(value_function, ComplexF64.(u), τ, derivs, derivs_t)
end

"""
То же, что θ, но с тэта-характеристикой, заданной через `char_shtiff`.

`faster` работает только тогда, когда `p` уже имеет тип `RiemannMatrix`.
"""
# !!! function θ_shtiff(u, p; derivs=[], derivs_t=[], char=zeros(length(u), length(u)), char_shtiff=[], faster=USE_FASTER_THETA, ϵ=1.0e-12)
function θ_shtiff(u, p; derivs=[], derivs_t=[], char=zeros(length(u), length(u)), char_shtiff=[], faster=USE_FASTER_THETA, ϵ=THETA_EPS) #???
    char_matrix = shtiff_char_matrix(length(u); char=char, char_shtiff=char_shtiff)
    char1, char2 = e_columns(char_matrix)
    return θ(u, p; derivs=derivs, derivs_t=derivs_t, char=[char1, char2], faster=faster, ϵ=ϵ)
end

"""
Вычисляет производную тэта-функции по элементу матрицы периода `Ω[i,j]`.
"""
# !!! function dθ_dΩ(z, p, i::Int, j::Int; char=nothing, derivs=[], derivs_t=[], faster=USE_FASTER_THETA, ϵ=1.0e-12)
function dθ_dΩ(z, p, i::Int, j::Int; char=nothing, derivs=[], derivs_t=[], faster=USE_FASTER_THETA, ϵ=THETA_EPS) #???
    derivs_t_full = vcat(collect(derivs_t), [[i, j]])
    if char === nothing
        return θ(z, p; derivs=derivs, derivs_t=derivs_t_full, faster=faster, ϵ=ϵ)
    end
    return θ(z, p; derivs=derivs, derivs_t=derivs_t_full, char=char, faster=faster, ϵ=ϵ)
end

"""
То же, что `dθ_dΩ`, но для специальной характеристики `char_shtiff`.
"""
# !!! function dθ_dΩ_shtiff(z, p, i::Int, j::Int; char=nothing, char_shtiff=[], derivs=[], derivs_t=[], faster=USE_FASTER_THETA, ϵ=1.0e-12)
function dθ_dΩ_shtiff(z, p, i::Int, j::Int; char=nothing, char_shtiff=[], derivs=[], derivs_t=[], faster=USE_FASTER_THETA, ϵ=THETA_EPS) #???
    char_matrix = if char === nothing
        shtiff_char_matrix(length(z); char_shtiff=char_shtiff)
    else
        shtiff_char_matrix(length(z); char=char, char_shtiff=char_shtiff)
    end
    char1, char2 = e_columns(char_matrix)
    return dθ_dΩ(z, p, i, j; char=[char1, char2], derivs=derivs, derivs_t=derivs_t, faster=faster, ϵ=ϵ)
end
