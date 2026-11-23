#??? Независимая проверка сокращения показателей без доступа к пакету Theta.jl.
using LinearAlgebra #???

module Theta #???
    struct RiemannMatrix #???
        τ::Matrix{ComplexF64} #???
    end #???
    RiemannMatrix(τ; siegel=false, ϵ=1e-12, nderivs=0) = #???
        RiemannMatrix(Matrix{ComplexF64}(τ)) #???
    function theta(z, R::RiemannMatrix; char=[]) #???
        a = isempty(char) ? zeros(2) : Float64.(char[1]) ./ 2 #???
        b = isempty(char) ? zeros(2) : Float64.(char[2]) ./ 2 #???
        result = 0.0 + 0.0im #???
        for k1 in -15:15, k2 in -15:15 #???
            n = [k1, k2] .+ a #???
            result += exp(π*im*(sum(n .* (R.τ*n)) + 2sum(n .* (z .+ b)))) #???
        end #???
        return result #???
    end #???
end #???

const THETA_EPS = 5e-16 #???
e_columns(matrix) = (matrix[:, 1], matrix[:, 2]) #???
special_chars() = ([0 0; 0 0], [1 0; 0 0], [1 1; 0 0], #???
                   [0 1; 1 0], [0 1; 1 1], [0 1; 0 1]) #???
function shtiff_char_matrix(g; char_shtiff=[]) #???
    result = zeros(Int, g, g) #???
    for i in char_shtiff #???
        result .+= special_chars()[i] #???
    end #???
    return result #???
end #???

include("Sigel.jl") #???
include("surface_theta_log.jl") #???
#??? Загружаем только изменённую функцию из func.jl: остальные пакеты здесь не нужны.
surface_definition = only(filter(e -> e isa Expr && e.head == :function && #???
    occursin("x_on_surface_new(", string(e.args[1])), #???
    Meta.parseall(read("func.jl", String)).args)) #???
eval(surface_definition) #???

function direct_surface_x(u, τ; zero_index=2, one_index=3, #???
                          infinity_index=6, added_index=1, sec_added_index=1) #???
    function θs(z, indices) #???
        char = shtiff_char_matrix(2; char_shtiff=indices) #???
        a, b = e_columns(char) #???
        return Theta.theta(z, Theta.RiemannMatrix(τ); char=[a, b]) #???
    end #???
    th1 = θs(u, [zero_index, added_index, 3, 5, sec_added_index]) #???
    th2 = θs(u, [infinity_index, added_index, 3, 5, sec_added_index]) #???
    th3 = θs([0, 0], [infinity_index, added_index, 3, 5, one_index]) #???
    th4 = θs([0, 0], [zero_index, added_index, 3, 5, one_index]) #???
    return -(th3/th4)^2 * (th1/th2)^2 #???
end #???

for (x, y, m, u, blocks) in ((0.4, 1.6, 0.3, ComplexF64[0.13+0.07im, 0.21+0.09im], :partial), #???
                              (0.4, 0.5, 0.3, ComplexF64[0.13+0.07im, 0.21+0.09im], :full)) #???
    τ = im .* ComplexF64[x m*x; m*x m*x+y] #???
    A, B, C, D = blocks == :partial ? #???
        ([0 0; 0 1], [1 0; 0 0], [-1 0; 0 0], [0 0; 0 1]) : #???
        (zeros(Int, 2, 2), Matrix{Int}(I, 2, 2), -Matrix{Int}(I, 2, 2), zeros(Int, 2, 2)) #???
    data = get_transform_data(τ, ComplexF64.(A), ComplexF64.(B), #???
                              ComplexF64.(C), ComplexF64.(D)) #???
    for indices in ([2, 1, 3, 5, 1], [6, 1, 3, 5, 1], #???
                    [6, 1, 3, 5, 3], [2, 1, 3, 5, 3]) #???
        via_theta = _surface_transformed_theta_log(u, data, indices) #???
        via_scaled_sum = _surface_transformed_theta_log(u, data, indices; use_theta=false) #???
        @assert abs(exp(via_theta - via_scaled_sum) - 1) < 1e-10 #???
    end #???
    for sec_added_index in (1, 2, 4, 5) #???
        actual = _surface_x_with_cancelled_exponents(u, τ, A, B, C, D; #???
            zero_index=2, one_index=3, infinity_index=6, #???
            added_index=1, sec_added_index=sec_added_index) #???
        expected = direct_surface_x(u, τ; sec_added_index=sec_added_index) #???
        error = abs(actual - expected) / max(1.0, abs(expected)) #???
        println(blocks, " sec=", sec_added_index, " error=", error) #???
        @assert error < 1e-8 #???
    end #???
end #???

for m in (1/11, 5/11), sec_added_index in (1, 2, 3, 4, 5, 6) #???
    τ = im .* ComplexF64[1e-5 m*1e-5; m*1e-5 m*1e-5+1e-5] #???
    A, B, C, D = choose_modular_blocks(τ) #???
    value = _surface_x_with_cancelled_exponents(ComplexF64[0, 0], τ, A, B, C, D; #???
        zero_index=2, one_index=3, infinity_index=6, #???
        added_index=5, sec_added_index=sec_added_index) #???
    println("small m=", m, " sec=", sec_added_index, " value=", value) #???
    @assert isfinite(real(value)) && isfinite(imag(value)) #???
    integrated_value = x_on_surface_new(ComplexF64[0, 0], τ; #???
        added_index=5, sec_added_index=sec_added_index) #???
    @assert isapprox(integrated_value, value; rtol=1e-14, atol=1e-14) #???
end #???
