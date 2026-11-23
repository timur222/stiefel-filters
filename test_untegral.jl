include("integral.jl")

N = 10
M = 6

x_bounds = (0.0001, 10.0)
y_bounds = (0.0001, 10.0)

x_steps = 40
y_steps = 40

proj_tol = 1e-8
proj_maxiter = 40
atol = 1e-5
rtol = 1e-4
absolute_jacobian = false
show_progress = true
progress_interval = 20

w(x, y, z) = exp(-x^2 - y^2 - z^2)

"""
Запускает эксперимент с расширенным интегралом и записывает результат в отдельный файл.
"""
function run_test_untegral(; output_path=joinpath(@__DIR__, "test_untegral_results.txt"))
    time_start = time()

    result = run_integral_case_extended(
        "extended_experiment";
        w=w,
        N_value=N,
        M_value=M,
        x_bounds=x_bounds,
        y_bounds=y_bounds,
        x_steps=x_steps,
        y_steps=y_steps,
        absolute_jacobian=absolute_jacobian,
        proj_tol=proj_tol,
        proj_maxiter=proj_maxiter,
        atol=atol,
        rtol=rtol,
        show_progress=show_progress,
        progress_interval=progress_interval,
    )

    open(output_path, "w") do io
        println(io, "Experiment: extended integral")
        println(io, "N = ", N)
        println(io, "M = ", M)
        println(io, "x_bounds = ", x_bounds)
        println(io, "y_bounds = ", y_bounds)
        println(io, "x_steps = ", x_steps)
        println(io, "y_steps = ", y_steps)
        println(io, "proj_tol = ", proj_tol)
        println(io, "proj_maxiter = ", proj_maxiter)
        println(io, "atol = ", atol)
        println(io, "rtol = ", rtol)
        println(io, "absolute_jacobian = ", absolute_jacobian)
        println(io, "show_progress = ", show_progress)
        println(io, "progress_interval = ", progress_interval)
        if result.error === nothing
            println(io, "integral = ", result.integral)
        else
            println(io, "error = ", result.error)
        end
        println(io, "elapsed = ", time() - time_start)
    end

    println("Results written to ", output_path)
    return result
end

script_path = joinpath(@__DIR__, basename(@__FILE__))
if !isempty(PROGRAM_FILE) && abspath(PROGRAM_FILE) == abspath(script_path)
    run_test_untegral()
end
