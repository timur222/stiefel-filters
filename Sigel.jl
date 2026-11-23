const MODULAR_TRANSFORM_CACHE = Dict{Any, Any}()
const MODULAR_CHAR_CACHE = Dict{Any, Any}()
const THETA_MATRIX_CACHE = Dict{Any, Theta.RiemannMatrix}() #???
const GENUS_ONE_NORMALIZER_CACHE = Dict{Float64, ComplexF64}() #???
#??? В сетке периода почти всегда различны; не храним все матрицы до конца построения.
const THETA_MATRIX_CACHE_LIMIT = 256 #???
const MODULAR_TRANSFORM_CACHE_LIMIT = 256 #???
const MODULAR_CHAR_CACHE_LIMIT = 4096 #???
const SIEGEL_TRANSFORM_THRESHOLD = 1e-3 #??? Применяем преобразование при диагональном элементе в 1000 раз меньше прежнего порога.

function _siegel_bounded_cache_store!(cache, key, value, limit) #???
    if length(cache) >= limit && !haskey(cache, key) #???
        empty!(cache) #???
    end #???
    cache[key] = value #???
    return value #???
end #???

"""
Возвращает матрицу периода τ: либо из RiemannMatrix, либо как копию обычной матрицы.
"""
function period_matrix(p)
    # !!! if p isa RiemannMatrix
    if p isa Theta.RiemannMatrix #???
        return Matrix{ComplexF64}(p.τ)
    end
    return Matrix{ComplexF64}(p)
end

"""
Выбирает вещественную величину, по которой решается, нужно ли делать модульное преобразование.

Если диагональный элемент почти чисто мнимый, используется его мнимая часть.
Иначе используется вещественная часть.
"""
function branch_value(x; atol=1e-12)
    if abs(imag(x)) > atol
        return imag(x)
    end
    return real(x)
end

"""
Возвращает блоки A, B, C, D для заранее заданного преобразования по правилу из проекта.
"""
function choose_modular_blocks(τ)
    g = size(τ, 1)
    A = Matrix{Int}(I, g, g)
    B = zeros(Int, g, g)
    C = zeros(Int, g, g)
    D = Matrix{Int}(I, g, g)

    if g != 2
        return A, B, C, D
    end

    p11 = branch_value(τ[1, 1])
    p22 = branch_value(τ[2, 2])

    # !!! if p11 > 1
    if p11 > SIEGEL_TRANSFORM_THRESHOLD #???
        return A, B, C, D
    # !!! elseif p11 < 1 && p22 > 1
    elseif p11 < SIEGEL_TRANSFORM_THRESHOLD && p22 > SIEGEL_TRANSFORM_THRESHOLD #???
        # !!! return [0 0; 0 1], [-1 0; 0 0], [1 0; 0 0], [0 0; 0 1]
        return [0 0; 0 1], [1 0; 0 0], [-1 0; 0 0], [0 0; 0 1] #??? S⁻¹ в первой координате, как и в ветви полного обращения.
    else
        return [0 0; 0 0], [1 0; 0 1], [-1 0; 0 -1], [0 0; 0 0]
    end
end

"""
Проверяет, что набор блоков задаёт тождественное преобразование.
"""
function is_identity_blocks(A, B, C, D)
    g = size(A, 1)
    return A == Matrix{Int}(I, g, g) &&
           B == zeros(Int, g, g) &&
           C == zeros(Int, g, g) &&
           D == Matrix{Int}(I, g, g)
end

"""
Нормализует характеристику: пустую заменяет на нулевую, непустую приводит к целым массивам.
"""
function normalize_char(char, g)
    if length(char) == 0
        return [zeros(Int, g), zeros(Int, g)]
    end
    return [Int.(char[1]), Int.(char[2])]
end

"""
Строит ключ кеша для массива по его содержимому.
"""
function matrix_cache_key(A)
    return Tuple(vec(ComplexF64.(A)))
end

#??? Переиспользуем матрицу Theta.jl для повторных вычислений при том же τ.
# !!! function theta_matrix_cached(τ; ϵ=1.0e-12)
# !!! function theta_matrix_cached(τ; ϵ=THETA_EPS)
function theta_matrix_cached(τ; ϵ=THETA_EPS, nderivs::Int=0) #??? Строим эллипсоиды только для требуемых производных.
    nderivs >= 0 || throw(ArgumentError("nderivs должен быть неотрицательным")) #???
    # !!! key = (size(τ), matrix_cache_key(τ), ϵ)
    key = (size(τ), matrix_cache_key(τ), ϵ, nderivs) #??? Производные требуют другого RiemannMatrix.
    #= !!! Прежний кеш рос без ограничений при построении большой сетки.
    return get!(THETA_MATRIX_CACHE, key) do
        Theta.RiemannMatrix(Matrix{ComplexF64}(τ), siegel=false, ϵ=ϵ)
    end
    =#
    haskey(THETA_MATRIX_CACHE, key) && return THETA_MATRIX_CACHE[key] #???
    # !!! matrix = Theta.RiemannMatrix(Matrix{ComplexF64}(τ), siegel=false, ϵ=ϵ)
    matrix = Theta.RiemannMatrix(Matrix{ComplexF64}(τ), siegel=false, ϵ=ϵ, nderivs=nderivs) #???
    # !!! return bounded_cache_store!(THETA_MATRIX_CACHE, key, matrix, THETA_MATRIX_CACHE_LIMIT)
    return _siegel_bounded_cache_store!(THETA_MATRIX_CACHE, key, matrix, THETA_MATRIX_CACHE_LIMIT) #???
end

"""
Строит ключ кеша для тэта-характеристики.
"""
function char_cache_key(char)
    return (Tuple(Int.(char[1])), Tuple(Int.(char[2])))
end

"""
Приводит блок модульного преобразования к целочисленному виду.
"""
function integer_block(A; atol=1e-10)
    A_int = Matrix{Int}(undef, size(A)...)
    for idx in eachindex(A)
        value = A[idx]
        if abs(imag(value)) > atol || abs(real(value) - round(real(value))) > atol
            error("Преобразование характеристик определено только для целочисленных блоков A, B, C, D.")
        end
        A_int[idx] = round(Int, real(value))
    end
    return A_int
end

"""
Приводит блоки A, B, C, D к целочисленному виду для работы с характеристиками.
"""
function integer_blocks(A, B, C, D; atol=1e-10)
    return (
        integer_block(A; atol=atol),
        integer_block(B; atol=atol),
        integer_block(C; atol=atol),
        integer_block(D; atol=atol),
    )
end

"""
Преобразует тэта-характеристику по блокам A, B, C, D.
"""
function transform_char(char, A, B, C, D)
    a = char[1]
    b = char[2]
    a_new = mod.(D * a - C * b + diag(C * transpose(D)), 2)
    b_new = mod.(-B * a + A * b + diag(A * transpose(B)), 2)
    return [Int.(a_new), Int.(b_new)]
end

"""
Вычисляет фазовый множитель, возникающий при преобразовании характеристики.
"""
function char_phase(char, A, B, C, D)
    if length(char) == 0
        return 1.0 + 0.0im
    end
    a = char[1]
    b = char[2]
    exponent = (transpose(a) * transpose(D) * B * a +
                transpose(b) * transpose(C) * A * b +
                2 * transpose(a) * transpose(B) * C * b)[1]
    return exp(π * im * exponent / 4)
end

"""
Применяет модульное преобразование к матрице периода τ.
"""
function transform_tau(τ, A, B, C, D)
    return (A * τ + B) * inv(C * τ + D)
end

"""
Редуцирует характеристику по модулю 2 и возвращает сопутствующий множитель.
"""
function reduce_char_with_factor(char)
    a = Int.(char[1])
    b = Int.(char[2])
    a_red = mod.(a, 2)
    b_red = mod.(b, 2)
    b_shift = div.(b - b_red, 2)
    factor = exp(π * im * sum(a_red .* b_shift))
    return [a_red, b_red], factor
end

"""
Подготавливает и кеширует данные преобразования, не зависящие от точки u и характеристики.
"""
function get_transform_data(τ, A, B, C, D)
    g = size(τ, 1)
    key = (
        g,
        matrix_cache_key(τ),
        matrix_cache_key(A),
        matrix_cache_key(B),
        matrix_cache_key(C),
        matrix_cache_key(D),
    )

    if haskey(MODULAR_TRANSFORM_CACHE, key)
        return MODULAR_TRANSFORM_CACHE[key]
    end

    M = C * τ + D
    M_inv = inv(M)
    τ_trans = transform_tau(τ, A, B, C, D)
    A_int, B_int, C_int, D_int = integer_blocks(A, B, C, D)
    Y_trans = Symmetric(0.5 .* (imag.(τ_trans) .+ transpose(imag.(τ_trans))))

    data = (
        transform_key=key,
        M_inv_T=transpose(M_inv),
        quad_matrix=M_inv * C,
        # !!! det_factor=det(M / im)^(-1 / 2),
        det_factor=det(D + im * C * τ)^(-1 / 2), #??? Правильная фазовая ветвь для частичного S⁻¹ и полного обращения.
        τ_trans=τ_trans,
        τ_trans_Y_inv=inv(Y_trans),
        τ_trans_lambda_min=minimum(eigvals(Y_trans)),
        A_int=A_int,
        B_int=B_int,
        C_int=C_int,
        D_int=D_int,
    )

    # !!! MODULAR_TRANSFORM_CACHE[key] = data
    # !!! bounded_cache_store!(MODULAR_TRANSFORM_CACHE, key, data, MODULAR_TRANSFORM_CACHE_LIMIT)
    _siegel_bounded_cache_store!(MODULAR_TRANSFORM_CACHE, key, data, MODULAR_TRANSFORM_CACHE_LIMIT) #???
    return data
end

"""
Подготавливает и кеширует данные преобразования, зависящие от характеристики.
"""
function get_char_data(transform_data, char_norm)
    char_red, char_factor = reduce_char_with_factor(char_norm)
    key = (transform_data.transform_key, char_cache_key(char_red), char_factor)

    if haskey(MODULAR_CHAR_CACHE, key)
        return MODULAR_CHAR_CACHE[key]
    end

    char_trans = transform_char(
        char_red,
        transform_data.A_int,
        transform_data.B_int,
        transform_data.C_int,
        transform_data.D_int,
    )
    g = length(char_trans[1])
    shift_n = 0.5 .* Float64.(char_trans[1])
    shift_x = 0.5 .* Float64.(char_trans[2])

    char_data = (
        char_trans=char_trans,
        shift_n=shift_n,
        shift_x=shift_x,
        g=g,
        base_factor=char_factor * char_phase(
            char_red,
            transform_data.A_int,
            transform_data.B_int,
            transform_data.C_int,
            transform_data.D_int,
        ) * transform_data.det_factor,
    )

    # !!! MODULAR_CHAR_CACHE[key] = char_data
    # !!! bounded_cache_store!(MODULAR_CHAR_CACHE, key, char_data, MODULAR_CHAR_CACHE_LIMIT)
    _siegel_bounded_cache_store!(MODULAR_CHAR_CACHE, key, char_data, MODULAR_CHAR_CACHE_LIMIT) #???
    return char_data
end

"""
Вычисляет тэта-функцию по прямой сумме по решётке для уже подготовленных данных.
"""
function theta_direct_reduced(z, transform_data, char_data; ϵ=1.0e-14)
    g = char_data.g
    y = Float64.(imag.(z))
    center = -(transform_data.τ_trans_Y_inv * y .+ char_data.shift_n)
    center_round = round.(Int, center)
    radius = max(
        8,
        ceil(Int, sqrt(max(0.0, -log(ϵ) / (π * transform_data.τ_trans_lambda_min))) + maximum(abs.(center .- center_round)) + 2),
    )

    total = 0.0 + 0.0im

    if g == 2
        τ11 = transform_data.τ_trans[1, 1]
        τ12 = transform_data.τ_trans[1, 2]
        τ22 = transform_data.τ_trans[2, 2]
        z1 = z[1] + char_data.shift_x[1]
        z2 = z[2] + char_data.shift_x[2]
        s1 = char_data.shift_n[1]
        s2 = char_data.shift_n[2]
        c1 = center_round[1]
        c2 = center_round[2]

        @inbounds for n1_int in (c1 - radius):(c1 + radius)
            n1 = n1_int + s1
            for n2_int in (c2 - radius):(c2 + radius)
                n2 = n2_int + s2
                quad = τ11 * (n1 * n1) + 2 * τ12 * n1 * n2 + τ22 * (n2 * n2)
                lin = 2 * (n1 * z1 + n2 * z2)
                total += exp(π * im * (quad + lin))
            end
        end
        return total
    end

    ranges = ntuple(i -> (center_round[i] - radius):(center_round[i] + radius), g)
    for idx in Iterators.product(ranges...)
        n = ComplexF64.(collect(idx)) .+ char_data.shift_n
        exponent = π * im * ((transpose(n) * transform_data.τ_trans * n)[1] + 2 * (transpose(n) * (z .+ char_data.shift_x))[1])
        total += exp(exponent)
    end

    return total
end

"""
Вычисляет значение тэта-функции после выбранного модульного преобразования.
"""
# !!! function transformed_theta_value(u, τ, A, B, C, D; char=[], ϵ=1.0e-14)
function transformed_theta_value(u, τ, A, B, C, D; char=[], ϵ=THETA_EPS) #???
    g = length(u)
    char_norm = normalize_char(char, g)
    transform_data = get_transform_data(τ, A, B, C, D)
    char_data = get_char_data(transform_data, char_norm)
    u_complex = ComplexF64.(u)
    z_trans = transform_data.M_inv_T * u_complex
    exponent = (transpose(u_complex) * transform_data.quad_matrix * u_complex)[1]
    prefactor = char_data.base_factor * exp(-π * im * exponent)
    # !!! return prefactor * theta_direct_reduced(ComplexF64.(vec(z_trans)), transform_data, char_data; ϵ=ϵ)
    return prefactor * Theta.theta(ComplexF64.(vec(z_trans)), theta_matrix_cached(transform_data.τ_trans; ϵ=ϵ); char=char_data.char_trans) #???
end

#= !!! Прежний вариант прямого суммирования без Theta.jl:
function theta_value_direct(u, τ; char=[], ϵ=1.0e-14)
    g = length(u)
    identity_block = Matrix{Int}(I, g, g)
    zero_block = zeros(Int, g, g)
    transform_data = get_transform_data(τ, identity_block, zero_block, zero_block, identity_block)
    char_red, char_factor = reduce_char_with_factor(normalize_char(char, g))
    char_data = (
        g=g,
        shift_n=0.5 .* Float64.(char_red[1]),
        shift_x=0.5 .* Float64.(char_red[2]),
    )
    return char_factor * theta_direct_reduced(ComplexF64.(u), transform_data, char_data; ϵ=ϵ)
end
=#

"""
Вычисляет значение тэта-функции без производных, автоматически применяя нужное модульное преобразование.
"""
# !!! function theta_value_auto(u, p; char=[], ϵ=1.0e-20)
function theta_value_auto(u, p; char=[], ϵ=THETA_EPS) #???
    τ = period_matrix(p)
    #??? В Theta.jl точки решётки рода 1 становятся скалярами, а затем складываются с вектором.
    #??? Для τ₂ = diag(τ, im) функция Theta.jl даёт θ(u|τ) * θ(0|im).
    if size(τ) == (1, 1) #???
        τ2 = ComplexF64[τ[1, 1] 0; 0 im] #???
        u2 = ComplexF64[u[1], 0] #???
        char2 = isempty(char) ? [] : [Int[char[1][1], 0], Int[char[2][1], 0]] #???
        normalizer = get!(GENUS_ONE_NORMALIZER_CACHE, Float64(ϵ)) do #???
            τ0 = ComplexF64[im 0; 0 im] #???
            sqrt(Theta.theta(ComplexF64[0, 0], theta_matrix_cached(τ0; ϵ=ϵ))) #???
        end #???
        return Theta.theta(u2, theta_matrix_cached(τ2; ϵ=ϵ); char=char2) / normalizer #???
    end #???
    A, B, C, D = choose_modular_blocks(τ)
    if is_identity_blocks(A, B, C, D)
        # !!! return theta(u, RiemannMatrix(τ, siegel=false, ϵ=ϵ), char=char)
        # !!! return theta_value_direct(u, τ; char=char, ϵ=ϵ)
        return Theta.theta(ComplexF64.(u), theta_matrix_cached(τ; ϵ=ϵ); char=char) #???
    end
    return transformed_theta_value(u, τ, ComplexF64.(A), ComplexF64.(B), ComplexF64.(C), ComplexF64.(D); char=char, ϵ=ϵ)
end
