#??? Вычисление x_on_surface_new после преобразования Зигеля без больших промежуточных экспонент.
#??? Общие множители det(Cτ+D)^(-1/2) и exp(-πiuᵀ(Cτ+D)^(-1)Cu)
#??? сокращаются в отношениях тета-функций до вычисления чисел.

const SURFACE_LOG_CONSTANT_CACHE = Ref{Any}(nothing) #???
const SURFACE_LOG_CONSTANT_LOCK = ReentrantLock() #???

function _surface_transformed_theta_log(u, transform_data, char_shtiff; #???
                                        ϵ=THETA_EPS, use_theta=true) #???
    g = length(u) #???
    g == 2 || throw(ArgumentError("логарифмическое вычисление поверхности поддерживает род 2")) #???
    char_matrix = shtiff_char_matrix(g; char_shtiff=char_shtiff) #???
    a, b = e_columns(char_matrix) #???
    char_red, char_factor = reduce_char_with_factor(normalize_char([a, b], g)) #???
    char_data = get_char_data(transform_data, [a, b]) #???
    modular_phase = char_factor * char_phase(char_red, transform_data.A_int, #???
                                            transform_data.B_int, transform_data.C_int, #???
                                            transform_data.D_int) #???

    z = transform_data.M_inv_T * ComplexF64.(u) #???
    τ = transform_data.τ_trans #???
    #??? В безопасном диапазоне оставляем вычисление ряда библиотеке Theta.jl.
    imag_z = imag.(z) #???
    theta_exp_log = π * dot(imag_z, transform_data.τ_trans_Y_inv * imag_z) #???
    if use_theta && transform_data.τ_trans_lambda_min > 0.01 && #???
       maximum(abs, imag.(τ)) < 100 && theta_exp_log < 300 #???
        value = Theta.theta(ComplexF64.(z), theta_matrix_cached(τ; ϵ=ϵ); #???
                            char=char_data.char_trans) #???
        if isfinite(real(value)) && isfinite(imag(value)) && !iszero(value) #???
            return log(modular_phase) + log(value) #???
        end #???
    end #???
    center = -(transform_data.τ_trans_Y_inv * imag.(z) .+ char_data.shift_n) #???
    nearest = round.(Int, center) #???
    λmin = transform_data.τ_trans_lambda_min #???
    λmin > 0 || throw(ArgumentError("мнимая часть преобразованной матрицы периода должна быть положительной")) #???
    #??? Из эллипсоида (n-center)ᵀY(n-center) получаем отдельный радиус для каждой оси.
    cutoff = -log(ϵ) / π #???
    radii = ntuple(j -> max(2, ceil(Int, sqrt(cutoff * #???
        transform_data.τ_trans_Y_inv[j, j]) + abs(center[j] - nearest[j]) + 2)), 2) #???

    z1 = z[1] + char_data.shift_x[1] #???
    z2 = z[2] + char_data.shift_x[2] #???
    τ11, τcross, τ22 = τ[1, 1], τ[1, 2] + τ[2, 1], τ[2, 2] #???
    highest_log = -Inf #???
    scaled_sum = 0.0 + 0.0im #???
    @inbounds for k1 in (nearest[1] - radii[1]):(nearest[1] + radii[1]) #???
        n1 = k1 + char_data.shift_n[1] #???
        for k2 in (nearest[2] - radii[2]):(nearest[2] + radii[2]) #???
            n2 = k2 + char_data.shift_n[2] #???
            exponent = π * im * (τ11*n1*n1 + τcross*n1*n2 + τ22*n2*n2 + #???
                                 2n1*z1 + 2n2*z2) #???
            log_magnitude = real(exponent) #???
            phase = cis(imag(exponent)) #???
            if log_magnitude > highest_log #???
                scaled_sum = scaled_sum * exp(highest_log - log_magnitude) + phase #???
                highest_log = log_magnitude #???
            else #???
                scaled_sum += exp(log_magnitude - highest_log) * phase #???
            end #???
        end #???
    end #???
    return log(modular_phase) + highest_log + log(scaled_sum) #???
end #???

function _surface_constant_log_ratio(transform_data, zero_index, one_index, #???
                                     infinity_index, added_index) #???
    key = (transform_data.transform_key, zero_index, one_index, infinity_index, #???
           added_index, THETA_EPS) #???
    cached = lock(SURFACE_LOG_CONSTANT_LOCK) do #???
        entry = SURFACE_LOG_CONSTANT_CACHE[] #???
        entry !== nothing && entry.key == key ? entry.value : nothing #???
    end #???
    cached !== nothing && return cached #???

    zero_u = zeros(ComplexF64, 2) #???
    numerator = _surface_transformed_theta_log(zero_u, transform_data, #???
        [infinity_index, added_index, 3, 5, one_index]) #???
    denominator = _surface_transformed_theta_log(zero_u, transform_data, #???
        [zero_index, added_index, 3, 5, one_index]) #???
    value = numerator - denominator #???
    lock(SURFACE_LOG_CONSTANT_LOCK) do #???
        SURFACE_LOG_CONSTANT_CACHE[] = (; key, value) #???
    end #???
    return value #???
end #???

function _surface_x_with_cancelled_exponents(u, τ, A, B, C, D; #???
                                             zero_index, one_index, infinity_index, #???
                                             added_index, sec_added_index) #???
    transform_data = get_transform_data(τ, ComplexF64.(A), ComplexF64.(B), #???
                                        ComplexF64.(C), ComplexF64.(D)) #???
    constant_log = _surface_constant_log_ratio(transform_data, zero_index, #???
                                               one_index, infinity_index, added_index) #???
    numerator_log = _surface_transformed_theta_log(u, transform_data, #???
        [zero_index, added_index, 3, 5, sec_added_index]) #???
    denominator_log = _surface_transformed_theta_log(u, transform_data, #???
        [infinity_index, added_index, 3, 5, sec_added_index]) #???
    return -exp(2 * (constant_log + numerator_log - denominator_log)) #???
end #???
