include("points.jl")
using QuadGK

"""
Форматирует число секунд в строку вида `hh:mm:ss`.
"""
function format_eta_seconds(seconds)
    total_seconds = max(0, round(Int, seconds))
    hours = total_seconds ÷ 3600
    minutes = (total_seconds % 3600) ÷ 60
    secs = total_seconds % 60
    return string(lpad(hours, 2, '0'), ":", lpad(minutes, 2, '0'), ":", lpad(secs, 2, '0'))
end

"""
Вычисляет интегранд после аналитически дешёвого интегрирования по `s`.
"""
function boundary_integrand_xy(w, x, y; N=N, M=M, absolute_jacobian=false, x4p_seed=nothing, return_seed=false, proj_tol=1e-10, proj_maxiter=100, atol=1e-9, rtol=1e-7)
    y45, next_seed = x_4_prime(
        x,
        y;
        N=N,
        M=M,
        active_seed=x4p_seed,
        return_seed=true,
        proj_tol=proj_tol,
        proj_maxiter=proj_maxiter,
    )
    y4 = x_4(x, y; N=N, M=M)
    y5 = x_5(x, y; N=N, M=M)
    y6 = x_6(x, y; N=N, M=M)

    inner_s, _ = quadgk(
        s -> w((1 - s) * y45 + s * y4, y5, y6),
        0.0,
        1.0;
        atol=atol,
        rtol=rtol,
    )

    jac = (y4 - y45) * (
        dx_x_5(x, y; N=N, M=M) * dy_x_6(x, y; N=N, M=M) -
        dy_x_5(x, y; N=N, M=M) * dx_x_6(x, y; N=N, M=M)
    )
    if absolute_jacobian
        jac = abs(jac)
    end

    value = inner_s * jac
    if return_seed
        return value, next_seed
    end
    return value
end

"""
Вычисляет интеграл из последней формулы по заданным границам `x` и `y`.

Интегрирование выполняется в порядке `s -> y -> x`, чтобы внутренний интеграл по `s`
не обращался к тэта-функциям повторно.
"""
function integrate_boundary_volume(w, x_bounds, y_bounds; N=N, M=M, absolute_jacobian=false, proj_tol=1e-10, proj_maxiter=100, atol=1e-9, rtol=1e-7)
    x_left, x_right = x_bounds
    y_left, y_right = y_bounds

    value, error = quadgk(
        x -> quadgk(
            y -> boundary_integrand_xy(
                w,
                x,
                y;
                N=N,
                M=M,
                absolute_jacobian=absolute_jacobian,
                proj_tol=proj_tol,
                proj_maxiter=proj_maxiter,
                atol=atol,
                rtol=rtol,
            ),
            y_left,
            y_right;
            atol=atol,
            rtol=rtol,
        )[1],
        x_left,
        x_right;
        atol=atol,
        rtol=rtol,
    )

    return value, error
end

"""
Вычисляет тот же интеграл по равномерной сетке по `x` и `y`, передавая seed
для `x_4_prime` от соседних точек.

Такой режим полезен, когда хочется стабилизировать `project_to_image` через continuation.
Интегрирование по `s` по-прежнему выполняется внутренним `quadgk`.
"""
function integrate_boundary_volume_seeded(w, x_bounds, y_bounds; N=N, M=M, x_steps=40, y_steps=40, absolute_jacobian=false, proj_tol=1e-8, proj_maxiter=40, atol=1e-9, rtol=1e-7)
    x_left, x_right = x_bounds
    y_left, y_right = y_bounds

    xs = collect(range(x_left, x_right; length=x_steps + 1))
    ys = collect(range(y_left, y_right; length=y_steps + 1))
    values = Matrix{Float64}(undef, length(xs), length(ys))

    column_seeds = fill(nothing, length(ys))

    for ix in eachindex(xs)
        row_seed = nothing
        new_column_seeds = Vector{Union{Nothing, ComplexF64}}(undef, length(ys))

        for iy in eachindex(ys)
            seed = row_seed === nothing ? column_seeds[iy] : row_seed
            value, next_seed = boundary_integrand_xy(
                w,
                xs[ix],
                ys[iy];
                N=N,
                M=M,
                absolute_jacobian=absolute_jacobian,
                x4p_seed=seed,
                return_seed=true,
                proj_tol=proj_tol,
                proj_maxiter=proj_maxiter,
                atol=atol,
                rtol=rtol,
            )
            values[ix, iy] = real(value)
            row_seed = next_seed
            new_column_seeds[iy] = next_seed
        end

        column_seeds = new_column_seeds
    end

    total = 0.0
    dx = (x_right - x_left) / x_steps
    dy = (y_right - y_left) / y_steps

    for ix in eachindex(xs)
        wx = (ix == 1 || ix == length(xs)) ? 0.5 : 1.0
        for iy in eachindex(ys)
            wy = (iy == 1 || iy == length(ys)) ? 0.5 : 1.0
            total += wx * wy * values[ix, iy]
        end
    end

    return total * dx * dy
end

"""
Вычисляет интегранд для расширенной области:
- первая полоса идёт по первой координате от `x_4_prime` до `x_4`,
- вторая полоса идёт по второй координате от `x_5` до `x_5_prime`
  при фиксированных `x_4` и `x_6`.
"""
function boundary_integrand_xy_extended(w, x, y;
    N=N,
    M=M,
    absolute_jacobian=false,
    x4p_seed=nothing,
    x5p_seed=nothing,
    return_seeds=false,
    proj_tol=1e-10,
    proj_maxiter=100,
    atol=1e-9,
    rtol=1e-7)

    y45, next_x4_seed = x_4_prime(
        x,
        y;
        N=N,
        M=M,
        active_seed=x4p_seed,
        return_seed=true,
        proj_tol=proj_tol,
        proj_maxiter=proj_maxiter,
    )
    y4 = x_4(x, y; N=N, M=M)
    y5 = x_5(x, y; N=N, M=M)
    y5prime, next_x5_seed = x_5_prime(
        x,
        y;
        N=N,
        M=M,
        active_seed=x5p_seed,
        return_seed=true,
        proj_tol=proj_tol,
        proj_maxiter=proj_maxiter,
    )
    y6 = x_6(x, y; N=N, M=M)

    inner_45, _ = quadgk(
        s -> w((1 - s) * y45 + s * y4, y5, y6),
        0.0,
        1.0;
        atol=atol,
        rtol=rtol,
    )
    inner_55, _ = quadgk(
        s -> w(y4, (1 - s) * y5 + s * y5prime, y6),
        0.0,
        1.0;
        atol=atol,
        rtol=rtol,
    )

    jac_45 = (y4 - y45) * (
        dx_x_5(x, y; N=N, M=M) * dy_x_6(x, y; N=N, M=M) -
        dy_x_5(x, y; N=N, M=M) * dx_x_6(x, y; N=N, M=M)
    )
    jac_55 = (y5prime - y5) * (
        dy_x_4(x, y; N=N, M=M) * dx_x_6(x, y; N=N, M=M) -
        dx_x_4(x, y; N=N, M=M) * dy_x_6(x, y; N=N, M=M)
    )

    if absolute_jacobian
        jac_45 = abs(jac_45)
        jac_55 = abs(jac_55)
    end

    value = inner_45 * jac_45 + inner_55 * jac_55
    if return_seeds
        return value, next_x4_seed, next_x5_seed
    end
    return value
end

"""
Вычисляет интеграл по расширенной области, где добавляется вторая полоса
`x_5 .. x_5_prime`.
"""
function integrate_boundary_volume_extended(w, x_bounds, y_bounds;
    N=N,
    M=M,
    absolute_jacobian=false,
    proj_tol=1e-10,
    proj_maxiter=100,
    atol=1e-9,
    rtol=1e-7)

    x_left, x_right = x_bounds
    y_left, y_right = y_bounds

    value, error = quadgk(
        x -> quadgk(
            y -> boundary_integrand_xy_extended(
                w,
                x,
                y;
                N=N,
                M=M,
                absolute_jacobian=absolute_jacobian,
                proj_tol=proj_tol,
                proj_maxiter=proj_maxiter,
                atol=atol,
                rtol=rtol,
            ),
            y_left,
            y_right;
            atol=atol,
            rtol=rtol,
        )[1],
        x_left,
        x_right;
        atol=atol,
        rtol=rtol,
    )

    return value, error
end

"""
Вычисляет расширенный интеграл по равномерной сетке по `x` и `y`, передавая
continuation seed отдельно для `x_4_prime` и `x_5_prime`.
"""
function integrate_boundary_volume_seeded_extended(w, x_bounds, y_bounds;
    N=N,
    M=M,
    x_steps=40,
    y_steps=40,
    absolute_jacobian=false,
    proj_tol=1e-8,
    proj_maxiter=40,
    atol=1e-9,
    rtol=1e-7,
    show_progress=false,
    progress_interval=10)

    x_left, x_right = x_bounds
    y_left, y_right = y_bounds

    xs = collect(range(x_left, x_right; length=x_steps + 1))
    ys = collect(range(y_left, y_right; length=y_steps + 1))
    values = Matrix{Float64}(undef, length(xs), length(ys))
    total_points = length(xs) * length(ys)
    processed_points = 0
    time_start = time()

    column_x4_seeds = fill(nothing, length(ys))
    column_x5_seeds = fill(nothing, length(ys))

    for ix in eachindex(xs)
        row_x4_seed = nothing
        row_x5_seed = nothing
        new_column_x4_seeds = Vector{Union{Nothing, ComplexF64}}(undef, length(ys))
        new_column_x5_seeds = Vector{Union{Nothing, ComplexF64}}(undef, length(ys))

        for iy in eachindex(ys)
            seed_x4 = row_x4_seed === nothing ? column_x4_seeds[iy] : row_x4_seed
            seed_x5 = row_x5_seed === nothing ? column_x5_seeds[iy] : row_x5_seed
            value, next_x4_seed, next_x5_seed = boundary_integrand_xy_extended(
                w,
                xs[ix],
                ys[iy];
                N=N,
                M=M,
                absolute_jacobian=absolute_jacobian,
                x4p_seed=seed_x4,
                x5p_seed=seed_x5,
                return_seeds=true,
                proj_tol=proj_tol,
                proj_maxiter=proj_maxiter,
                atol=atol,
                rtol=rtol,
            )
            values[ix, iy] = real(value)
            row_x4_seed = next_x4_seed
            row_x5_seed = next_x5_seed
            new_column_x4_seeds[iy] = next_x4_seed
            new_column_x5_seeds[iy] = next_x5_seed

            processed_points += 1
            if show_progress && ((processed_points % progress_interval == 0) || (processed_points == total_points))
                elapsed = time() - time_start
                remaining = elapsed * (total_points - processed_points) / processed_points
                percent = round(100 * processed_points / total_points; digits=1)
                println(
                    "progress ",
                    processed_points,
                    "/",
                    total_points,
                    " (",
                    percent,
                    "%), elapsed=",
                    format_eta_seconds(elapsed),
                    ", eta=",
                    format_eta_seconds(remaining),
                )
            end
        end

        column_x4_seeds = new_column_x4_seeds
        column_x5_seeds = new_column_x5_seeds
    end

    total = 0.0
    dx = (x_right - x_left) / x_steps
    dy = (y_right - y_left) / y_steps

    for ix in eachindex(xs)
        wx = (ix == 1 || ix == length(xs)) ? 0.5 : 1.0
        for iy in eachindex(ys)
            wy = (iy == 1 || iy == length(ys)) ? 0.5 : 1.0
            total += wx * wy * values[ix, iy]
        end
    end

    return total * dx * dy
end

"""
Запускает один тест интеграла и возвращает именованный кортеж с результатом или ошибкой.
"""
function run_integral_case(name; N_value=N, M_value=M, x_bounds=(0.0001, 0.905), y_bounds=(0.0001, 0.705), atol=1e-5, rtol=1e-4)
    time_int = time()
    try
        integral = integrate_boundary_volume_seeded(
            (x, y, z) -> exp(-x^2 - y^2 - z^2),
            x_bounds,
            y_bounds;
            N=N_value,
            M=M_value,
            atol=atol,
            rtol=rtol,
        )
        return (
            name=name,
            N=N_value,
            M=M_value,
            x_bounds=x_bounds,
            y_bounds=y_bounds,
            integral=integral,
            elapsed=time() - time_int,
            error=nothing,
        )
    catch err
        return (
            name=name,
            N=N_value,
            M=M_value,
            x_bounds=x_bounds,
            y_bounds=y_bounds,
            integral=nothing,
            elapsed=time() - time_int,
            error=sprint(showerror, err),
        )
    end
end

"""
Запускает один тест расширенного интеграла и возвращает именованный кортеж с результатом или ошибкой.
"""
function run_integral_case_extended(name;
    w=(x, y, z) -> exp(-x^2 - y^2 - z^2),
    N_value=N,
    M_value=M,
    x_bounds=(0.0001, 0.905),
    y_bounds=(0.0001, 0.705),
    x_steps=40,
    y_steps=40,
    absolute_jacobian=false,
    proj_tol=1e-8,
    proj_maxiter=40,
    atol=1e-5,
    rtol=1e-4,
    show_progress=false,
    progress_interval=10)

    time_int = time()
    try
        integral = integrate_boundary_volume_seeded_extended(
            w,
            x_bounds,
            y_bounds;
            N=N_value,
            M=M_value,
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
        return (
            name=name,
            N=N_value,
            M=M_value,
            x_bounds=x_bounds,
            y_bounds=y_bounds,
            x_steps=x_steps,
            y_steps=y_steps,
            absolute_jacobian=absolute_jacobian,
            proj_tol=proj_tol,
            proj_maxiter=proj_maxiter,
            integral=integral,
            elapsed=time() - time_int,
            error=nothing,
        )
    catch err
        return (
            name=name,
            N=N_value,
            M=M_value,
            x_bounds=x_bounds,
            y_bounds=y_bounds,
            x_steps=x_steps,
            y_steps=y_steps,
            absolute_jacobian=absolute_jacobian,
            proj_tol=proj_tol,
            proj_maxiter=proj_maxiter,
            integral=nothing,
            elapsed=time() - time_int,
            error=sprint(showerror, err),
        )
    end
end

"""
Записывает результаты тестов интеграла в текстовый файл.
"""
function write_integral_results(filepath; cases=[])
    open(filepath, "w") do io
        for case in cases
            println(io, "Case: ", case.name)
            println(io, "N = ", case.N, ", M = ", case.M)
            println(io, "x_bounds = ", case.x_bounds)
            println(io, "y_bounds = ", case.y_bounds)
            if case.error === nothing
                println(io, "Integral = ", case.integral)
            else
                println(io, "Error = ", case.error)
            end
            println(io, "Time = ", case.elapsed)
            println(io)
        end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    u = [0.1 + 0.7 * im, 0.2 + 0.8 * im]
    du = [0.1 + 0.2 * im, 0.3 + 0.4 * im]
    p = [0.5 + 1.2 * im 0.3 + 0.4 * im; 0.3 + 0.4 * im 0.6 + 1.5 * im]
    dp = [0.1 + 0.8 * im 0.3 + 0.4 * im; 0.3 + 0.4 * im 0.5 + 0.6 * im]

    θ(u, p)

    time_reg = time()
    theta_val = θ(du, dp)
    println("θ = ", norm(theta_val))
    println("time_reg = ", time() - time_reg)

    cases = [
        run_integral_case(
            "M4_small";
            N_value=10,
            M_value=4,
            x_bounds=(0.0001, 0.905),
            y_bounds=(0.0001, 0.705),
            atol=1e-5,
            rtol=1e-4,
        ),
        run_integral_case(
            "M4_large";
            N_value=10,
            M_value=4,
            x_bounds=(0.0001, 100),
            y_bounds=(0.0001, 100),
            atol=1e-5,
            rtol=1e-4,
        ),
        run_integral_case(
            "M6_large";
            N_value=10,
            M_value=6,
            x_bounds=(0.0001, 100),
            y_bounds=(0.0001, 100),
            atol=1e-5,
            rtol=1e-4,
        ),
    ]

    output_path = joinpath(@__DIR__, "integral_results.txt")
    #write_integral_results(output_path; cases=cases)
    #println("Integral results written to ", output_path)

    #test_der_x_on_surface_new(u, p, du, dp)
end

