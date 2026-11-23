#??? Независимая проверка модульного преобразования без установленного пакета Theta.jl.
using LinearAlgebra #???

module Theta #???
    struct RiemannMatrix #???
        τ::Matrix{ComplexF64} #???
        nderivs::Int #???
    end #???
    RiemannMatrix(τ; siegel=false, ϵ=1e-12, nderivs=0) = #???
        RiemannMatrix(Matrix{ComplexF64}(τ), nderivs) #???

    function theta(z, R::RiemannMatrix; char=[]) #???
        a = isempty(char) ? zeros(2) : Float64.(char[1]) ./ 2 #???
        b = isempty(char) ? zeros(2) : Float64.(char[2]) ./ 2 #???
        total = 0.0 + 0.0im #???
        for n1 in -12:12, n2 in -12:12 #???
            n = Float64[n1, n2] .+ a #???
            exponent = π * im * (sum(n .* (R.τ * n)) + 2sum(n .* (z .+ b))) #???
            total += exp(exponent) #???
        end #???
        return total #???
    end #???
end #???

const THETA_EPS = 5e-16 #???
include("Sigel.jl") #???

function check_case(name, τ, A, B, C, D) #???
    worst_error = 0.0 #???
    worst_char = nothing #???
    worst_ratio = NaN + NaN*im #???
    chars = [[Int[a1,a2], Int[b1,b2]] #???
             for a1 in 0:1 for a2 in 0:1 for b1 in 0:1 for b2 in 0:1] #???
    append!(chars, ([[2,1], [1,0]], [[3,-1], [2,1]], #???
                    [[-2,0], [1,-2]], [[1,3], [-1,2]])) #??? Характеристики из сумм особых точек.
    for z in (ComplexF64[0.13+0.07im, 0.21+0.09im], #???
              ComplexF64[0.5, 0.5], ComplexF64[0.5-0.2im, 0.5-0.1im]) #???
        # !!! for a1 in 0:1, a2 in 0:1, b1 in 0:1, b2 in 0:1
        for char in chars #???
            expected = Theta.theta(z, Theta.RiemannMatrix(τ); char=char) #???
            abs(expected) < 1e-9 && continue #???
            # !!! actual = transformed_theta_value(z, τ, ComplexF64.(A), ComplexF64.(B),
            # !!!                                  ComplexF64.(C), ComplexF64.(D); char=char)
            actual = theta_value_auto(z, τ; char=char) #??? Проверяем именно общедоступный путь вызова.
            relative_error = abs(actual-expected) / abs(expected) #???
            if relative_error > worst_error #???
                worst_error = relative_error #???
                worst_char = char #???
                worst_ratio = actual / expected #???
            end #???
        end #???
    end #???
    println(name, " worst_error=", worst_error, " char=", worst_char, #???
            " ratio=", worst_ratio) #???
    @assert worst_error < 1e-9 "Ошибка модульного преобразования в случае $name" #???
end #???

for (name, x, y, m) in (("identity", 1.4, 0.8, 0.3), #???
                        ("partial", 0.4, 1.6, 0.3), #???
                        ("full", 0.4, 0.4, 0.3), #???
                        ("partial near boundary", 0.999, 1.2, 1/11), #???
                        ("full near boundary", 0.4, 0.88, 1/11), #???
                        ("partial other M", 0.2, 1.3, 5/11)) #???
    τ = ComplexF64[x m*x; m*x m*x+y] .* im #???
    A, B, C, D = choose_modular_blocks(τ) #???
    symplectic_error = maximum(abs.(transpose(A)*D-transpose(C)*B-I)) #???
    symmetry_error = maximum(abs.(transform_tau(τ,A,B,C,D) - transpose(transform_tau(τ,A,B,C,D)))) #???
    minimum_imaginary_eigenvalue = minimum(eigvals(Symmetric(imag.(transform_tau(τ,A,B,C,D))))) #???
    println(name, " symplectic_error=", symplectic_error, #???
            " symmetry_error=", symmetry_error, #???
            " transformed_min_imag=", minimum_imaginary_eigenvalue) #???
    # !!! is_identity_blocks(A,B,C,D) || check_case(name, τ, A,B,C,D)
    check_case(name, τ, A,B,C,D) #???
end #???

for (name, τ) in (("partial complex", ComplexF64[0.1+0.4im 0.05+0.12im; #???
                                                  0.05+0.12im 0.2+1.7im]), #???
                  ("full complex", ComplexF64[0.1+0.4im 0.05+0.12im; #???
                                               0.05+0.12im 0.2+0.52im])) #???
    A,B,C,D = choose_modular_blocks(τ) #???
    check_case(name, τ, A,B,C,D) #???
end #???
