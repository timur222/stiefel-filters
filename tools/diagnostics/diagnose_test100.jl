# !!! Standalone diagnostic: runs the project's custom theta path with stdlib only.
#??? Diagnostic uses the installed Theta.jl package.
# !!! Run with JULIA_LOAD_PATH=@stdlib julia --startup-file=no --compiled-modules=no diagnose_test100.jl
#??? Run with julia --startup-file=no diagnose_test100.jl
module DiagnoseTest100 #??? Keep diagnostic definitions out of Main.
using LinearAlgebra
using Printf
# !!! abstract type RiemannMatrix end # Only used in isa checks; p is a plain matrix here.
using Theta #???
include("../../src/stiefel/theta.jl")
include("../../src/stiefel/surface_theta_log.jl")
core_file = joinpath(@__DIR__, "../../src/stiefel/func.jl")
source = first(split(read(core_file, String), "function shtiff_fraction(t,"))
source = replace(source, "using Polynomials" => "", "import Makie" => "",
                 "include(\"theta.jl\")" => "", "include(\"surface_theta_log.jl\")" => "")
include_string(@__MODULE__, source, core_file)

# Polynomial construction and Horner evaluation, without the plotting dependencies.
function fraction(z, p, c)
    function coefficients(roots)
        a = ComplexF64[1]
        for r in roots
            b = zeros(ComplexF64, length(a) + 1)
            b[1:end-1] .-= r .* a
            b[2:end] .+= a
            a = b
        end
        a
    end
    a, b = c .* coefficients(z), coefficients(p)
    x -> evalpoly(x, a) / evalpoly(x, b)
end

function run_case(epsilon)
    p = im .* [1.0 0.5; 0.5 0.5 * (1 + epsilon)]
    println("\nEPS=", epsilon)
    @show choose_modular_blocks(p)
    bottom = path_on_jacobian_new(0., 0., 0.5im, p)
    @show bottom θ_shtiff([-0.5im, bottom], p; char_shtiff=[3,5])
    result = shtiff_fraction_new(0.2, 10, p, 5)
    @show result
    z, poles, c = -result[1], -result[2], -result[3]
    expanded = fraction(z, poles, c)
    factored = x -> c * prod(x .- z) / prod(x .- poles)
    @show maximum(abs(expanded(x) - factored(x)) for x in range(-2., 2.; length=10001))
    @show maximum(abs(imag(factored(x))) for x in range(-2., 2.; length=10001))
    for k in 0:10
        target = 0.25 - 0.05im*k
        up = path_on_jacobian_new(0., 0., target, p)
        down = path_on_jacobian_new(-0.5im, bottom, target, p)
        @printf("k=%2d up=%+.10f%+.10fi down=%+.10f%+.10fi x_up=%s x_down=%s\n", k, real(up), imag(up), real(down), imag(down), string(x_on_surface_new([target,up],p)), string(x_on_surface_new([target,down],p)))
    end
end
for epsilon in (0.9, 0.8)
    run_case(epsilon)
end

# Change only the sign of the starting endpoint, in an in-memory variant.
fixed_source = source[findfirst("function shtiff_fraction_new", source).start:end]
fixed_source = replace(fixed_source,
    "function shtiff_fraction_new(" => "function shtiff_fraction_new_fixed(",
    "u_curr[2], 0.25*im*n*t, p)" => "u_curr[2], -0.25*im*n*t, p)")
# !!! include_string(Main, fixed_source, "sign_fix_experiment")
include_string(@__MODULE__, fixed_source, "sign_fix_experiment") #???
println("\nSIGN FIX: source files unchanged")
for epsilon in 0.8:0.01:0.9
    p = im .* [1.0 0.5; 0.5 0.5 * (1 + epsilon)]
    bottom = path_on_jacobian_new(0., 0., -0.5im, p)
    result = redirect_stdout(devnull) do
        shtiff_fraction_new_fixed(0.2, 10, p, 5)
    end
    @printf("eps=%.2f bottom=%s residual=%.3g coefficient=%s\n", epsilon, string(bottom), abs(θ_shtiff([-0.5im,bottom],p;char_shtiff=[3,5])), string(result[3]))
    if epsilon == 0.8
        @show result
    end
end

trace_source = source[findfirst("function path_on_jacobian_new", source).start:findfirst("function shtiff_fraction_new", source).start-1]
trace_source = replace(trace_source, "path_on_jacobian_new" => "trace_path",
    "u_21_curr = u_proj[1]" => "push!(path_trace, (start=ComplexF64[u_21_prev,u_22_prev], trial=u_trial, projected=u_proj, residual=residual)); u_21_curr = u_proj[1]")
# !!! include_string(Main, trace_source, "trace_experiment")
include_string(@__MODULE__, trace_source, "trace_experiment") #???
const path_trace = Any[]
for sign in (1,-1)
    p = im .* [1.0 0.5; 0.5 0.9]
    bottom = path_on_jacobian_new(0.,0.,sign*0.5im,p)
    empty!(path_trace)
    trace_path(-0.5im,bottom,0.25-0.3im,p)
    distances = [norm(v.projected-v.start) for v in path_trace]
    @show sign path_trace[1] maximum(distances) path_trace[argmax(distances)]
end
end #??? module DiagnoseTest100
