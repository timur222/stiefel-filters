include("points.jl")
using QuadGK

"""
Вычисляет производную `x_on_surface_new` по направлениям `du` и `dp`.

Если `p` - обычная матрица, внутри автоматически применяется то же модульное преобразование,
что и в базовой функции `θ`.
"""
function der_x_on_surface_new(u, p, du, dp; zero_index=3, one_index=2, infinity_index=1, added_index=4)
    s = length(u)
    use_faster = p isa RiemannMatrix

    function deriv_u_theta_shtiff(u_point, char_shtiff)
        return sum(
            du[i] * θ_shtiff(
                u_point,
                p;
                derivs=[[(j == i) ? 1 : 0 for j in 1:s]],
                char_shtiff=char_shtiff,
                faster=use_faster,
            )
            for i in 1:s
        )
    end

    function deriv_p_theta_shtiff(u_point, char_shtiff)
        total = zero(ComplexF64)
        for i in 1:s
            for j in i:s
                total += dp[i, j] * dθ_dΩ_shtiff(
                    u_point,
                    p,
                    i,
                    j;
                    char_shtiff=char_shtiff,
                    faster=use_faster,
                )
            end
        end
        return total
    end

    char1 = [zero_index, added_index, 3, 5]
    char2 = [infinity_index, added_index, 3, 5]
    char3 = [infinity_index, added_index, 3, 5, one_index]
    char4 = [zero_index, added_index, 3, 5, one_index]

    zero_u = zeros(ComplexF64, s)

    th1 = θ_shtiff(u, p; char_shtiff=char1, faster=use_faster)
    th2 = θ_shtiff(u, p; char_shtiff=char2, faster=use_faster)
    th3 = θ_shtiff(zero_u, p; char_shtiff=char3, faster=use_faster)
    th4 = θ_shtiff(zero_u, p; char_shtiff=char4, faster=use_faster)

    dth1 = deriv_u_theta_shtiff(u, char1) + deriv_p_theta_shtiff(u, char1)
    dth2 = deriv_u_theta_shtiff(u, char2) + deriv_p_theta_shtiff(u, char2)
    dth3 = deriv_p_theta_shtiff(zero_u, char3)
    dth4 = deriv_p_theta_shtiff(zero_u, char4)

    x_val = x_on_surface_new(
        u,
        p;
        zero_index=zero_index,
        one_index=one_index,
        infinity_index=infinity_index,
        added_index=added_index,
    )

    return 2 * x_val * (dth3 / th3 - dth4 / th4 + dth1 / th1 - dth2 / th2)
end

"""
Вычисляет производную проекции `project_to_image` по направлениям `du` и `dp`.

Активные координаты определяются из условий `theta = 0` через неявную функцию,
а неактивные координаты сохраняют производные исходной точки.
"""
function der_project_to_image(u, p, du, dp; shtiff_chars=[[]], act_vars=[1], tol=1e-10, maxiter=100)
    s = length(u)
    act_vars = collect(act_vars)
    shtiff_chars = collect(shtiff_chars)
    use_faster = p isa RiemannMatrix

    if isempty(shtiff_chars) || all(isempty, shtiff_chars)
        return copy(u), ComplexF64.(du)
    end

    u_proj = project_to_image(u, p; shtiff_chars=shtiff_chars, act_vars=act_vars, tol=tol, maxiter=maxiter)
    du_proj = zeros(ComplexF64, s)

    act_set = Set(act_vars)
    fixed_vars = [j for j in 1:s if !(j in act_set)]
    for j in fixed_vars
        du_proj[j] = du[j]
    end

    rhs = zeros(ComplexF64, length(shtiff_chars))
    for (char_index, char_k) in enumerate(shtiff_chars)
        fixed_part = sum(
            du[j] * θ_shtiff(
                u_proj,
                p;
                char_shtiff=char_k,
                derivs=[[(idx == j) ? 1 : 0 for idx in 1:s]],
                faster=use_faster,
            )
            for j in fixed_vars
        )

        period_part = sum(
            dp[i, j] * dθ_dΩ_shtiff(
                u_proj,
                p,
                i,
                j;
                char_shtiff=char_k,
                faster=use_faster,
            )
            for i in 1:s for j in i:s
        )

        rhs[char_index] = fixed_part + period_part
    end

    jac = Matrix{ComplexF64}(undef, length(shtiff_chars), length(act_vars))
    for (char_index, char_k) in enumerate(shtiff_chars)
        for (act_index, global_index) in enumerate(act_vars)
            jac[char_index, act_index] = θ_shtiff(
                u_proj,
                p;
                char_shtiff=char_k,
                derivs=[[(idx == global_index) ? 1 : 0 for idx in 1:s]],
                faster=use_faster,
            )
        end
    end

    act_der = -(jac \ rhs)
    for (local_index, global_index) in enumerate(act_vars)
        du_proj[global_index] = act_der[local_index]
    end

    return u_proj, du_proj
end

"""
Проверяет производную `der_x_on_surface_new` через разностное отношение.
"""
function test_der_x_on_surface_new(u, p, du, dp;
    zero_index=3, one_index=2, infinity_index=1, added_index=4,
    steps=(1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-11, 1e-12, 1e-13, 1e-14, 1e-15))

    analytic = der_x_on_surface_new(
        u,
        p,
        du,
        dp;
        zero_index=zero_index,
        one_index=one_index,
        infinity_index=infinity_index,
        added_index=added_index,
    )

    x0 = x_on_surface_new(
        u,
        p;
        zero_index=zero_index,
        one_index=one_index,
        infinity_index=infinity_index,
        added_index=added_index,
    )

    println("analytic = ", analytic)
    for s in steps
        diff_quot = (
            x_on_surface_new(
                u .+ s .* du,
                p .+ s .* dp;
                zero_index=zero_index,
                one_index=one_index,
                infinity_index=infinity_index,
                added_index=added_index,
            ) - x0
        ) / s
        println("s = ", s, " err = ", abs(diff_quot - analytic))
    end
end

"""
Строит матрицу периода `p=iΩ` и её производные по `x` и `y`.
"""
function xy_period_data(x, y; N=N, M=M)
    m = M / N
    OM = ComplexF64[x m * x; m * x m * x + y]
    p = OM * im
    dp_dx = ComplexF64[1 m; m m] * im
    dp_dy = ComplexF64[0 0; 0 1] * im
    return OM, p, dp_dx, dp_dy
end

"""
Возвращает данные для `x_4`.
"""
function x_4_data(x, y; N=N, M=M)
    OM, p, dp_dx, dp_dy = xy_period_data(x, y; N=N, M=M)
    m = M / N
    u = ComplexF64[0.5 - OM[1, 2] * im / 2, -im * OM[2, 2] / 2]
    du_dx = ComplexF64[-im * m / 2, -im * m / 2]
    du_dy = ComplexF64[0, -im / 2]
    return u, p, du_dx, du_dy, dp_dx, dp_dy
end

"""
Возвращает данные для `x_5`.
"""
function x_5_data(x, y; N=N, M=M)
    OM, p, dp_dx, dp_dy = xy_period_data(x, y; N=N, M=M)
    m = M / N
    u = ComplexF64[0.5 - OM[1, 2] * im / 2, 0.5 - im * OM[2, 2] / 2]
    du_dx = ComplexF64[-im * m / 2, -im * m / 2]
    du_dy = ComplexF64[0, -im / 2]
    return u, p, du_dx, du_dy, dp_dx, dp_dy
end

"""
Возвращает данные для `x_6`.
"""
function x_6_data(x, y; N=N, M=M)
    _, p, dp_dx, dp_dy = xy_period_data(x, y; N=N, M=M)
    u = ComplexF64[0.5, 0.5]
    du_dx = zeros(ComplexF64, 2)
    du_dy = zeros(ComplexF64, 2)
    return u, p, du_dx, du_dy, dp_dx, dp_dy
end

"""
Возвращает данные для `x_4_prime`, включая производную проекции на `theta[35]=0`.
"""
function x_4_prime_projected_data(x, y; N=N, M=M, tol=1e-10, maxiter=100)
    OM, p, dp_dx, dp_dy = xy_period_data(x, y; N=N, M=M)
    alpha = (M + 1) / M
    m = M / N

    u = ComplexF64[
        0.5 - OM[1, 2] * im / 2 * alpha,
        -im * OM[2, 2] / 2 * alpha,
    ]
    du_dx = ComplexF64[-im * m * alpha / 2, -im * m * alpha / 2]
    du_dy = ComplexF64[0, -im * alpha / 2]

    u_proj_dx, du_proj_dx = der_project_to_image(
        u,
        p,
        du_dx,
        dp_dx;
        shtiff_chars=[[3, 5]],
        act_vars=[2],
        tol=tol,
        maxiter=maxiter,
    )
    u_proj_dy, du_proj_dy = der_project_to_image(
        u,
        p,
        du_dy,
        dp_dy;
        shtiff_chars=[[3, 5]],
        act_vars=[2],
        tol=tol,
        maxiter=maxiter,
    )

    return u_proj_dx, du_proj_dx, du_proj_dy, p, dp_dx, dp_dy
end

"""
Вычисляет частную производную `x_4` по `x`.
"""
function dx_x_4(x, y; N=N, M=M)
    u, p, du_dx, _, dp_dx, _ = x_4_data(x, y; N=N, M=M)
    return real(-der_x_on_surface_new(u, p, du_dx, dp_dx; added_index=5))
end

"""
Вычисляет частную производную `x_4` по `y`.
"""
function dy_x_4(x, y; N=N, M=M)
    u, p, _, du_dy, _, dp_dy = x_4_data(x, y; N=N, M=M)
    return real(-der_x_on_surface_new(u, p, du_dy, dp_dy; added_index=5))
end

"""
Вычисляет частную производную `x_5` по `x`.
"""
function dx_x_5(x, y; N=N, M=M)
    u, p, du_dx, _, dp_dx, _ = x_5_data(x, y; N=N, M=M)
    return real(-der_x_on_surface_new(u, p, du_dx, dp_dx; added_index=6))
end

"""
Вычисляет частную производную `x_5` по `y`.
"""
function dy_x_5(x, y; N=N, M=M)
    u, p, _, du_dy, _, dp_dy = x_5_data(x, y; N=N, M=M)
    return real(-der_x_on_surface_new(u, p, du_dy, dp_dy; added_index=6))
end

"""
Вычисляет частную производную `x_6` по `x`.
"""
function dx_x_6(x, y; N=N, M=M)
    u, p, du_dx, _, dp_dx, _ = x_6_data(x, y; N=N, M=M)
    return real(-der_x_on_surface_new(u, p, du_dx, dp_dx; added_index=4))
end

"""
Вычисляет частную производную `x_6` по `y`.
"""
function dy_x_6(x, y; N=N, M=M)
    u, p, _, du_dy, _, dp_dy = x_6_data(x, y; N=N, M=M)
    return real(-der_x_on_surface_new(u, p, du_dy, dp_dy; added_index=4))
end

"""
Вычисляет частную производную `x_4_prime` по `x`.
"""
function dx_x_4_prime(x, y; N=N, M=M, tol=1e-10, maxiter=100)
    u_proj, du_proj_dx, _, p, dp_dx, _ = x_4_prime_projected_data(x, y; N=N, M=M, tol=tol, maxiter=maxiter)
    return real(-der_x_on_surface_new(u_proj, p, du_proj_dx, dp_dx; added_index=5))
end

"""
Вычисляет частную производную `x_4_prime` по `y`.
"""
function dy_x_4_prime(x, y; N=N, M=M, tol=1e-10, maxiter=100)
    u_proj, _, du_proj_dy, p, _, dp_dy = x_4_prime_projected_data(x, y; N=N, M=M, tol=tol, maxiter=maxiter)
    return real(-der_x_on_surface_new(u_proj, p, du_proj_dy, dp_dy; added_index=5))
end

"""
Проверяет частные производные функции `f(x,y)` по `x` и `y` через разностные отношения.
"""
function test_xy_derivative(f, fx, fy, x, y; steps=(1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7))
    analytic_x = fx(x, y)
    analytic_y = fy(x, y)
    f0 = f(x, y)

    println("analytic_x = ", analytic_x)
    for s in steps
        diff_quot = (f(x + s, y) - f0) / s
        println("s = ", s, " err_x = ", abs(diff_quot - analytic_x))
    end

    println("analytic_y = ", analytic_y)
    for s in steps
        diff_quot = (f(x, y + s) - f0) / s
        println("s = ", s, " err_y = ", abs(diff_quot - analytic_y))
    end
end

"""
Проверяет производные `x_4`, `x_4_prime`, `x_5`, `x_6` через разность.
"""
function test_boundary_derivatives(x, y; N=N, M=M, steps=(1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7))
    println("x_4")
    test_xy_derivative(
        (x_curr, y_curr) -> x_4(x_curr, y_curr; N=N, M=M),
        (x_curr, y_curr) -> dx_x_4(x_curr, y_curr; N=N, M=M),
        (x_curr, y_curr) -> dy_x_4(x_curr, y_curr; N=N, M=M),
        x,
        y;
        steps=steps,
    )
    println()

    println("x_4_prime")
    test_xy_derivative(
        (x_curr, y_curr) -> x_4_prime(x_curr, y_curr; N=N, M=M),
        (x_curr, y_curr) -> dx_x_4_prime(x_curr, y_curr; N=N, M=M),
        (x_curr, y_curr) -> dy_x_4_prime(x_curr, y_curr; N=N, M=M),
        x,
        y;
        steps=steps,
    )
    println()

    println("x_5")
    test_xy_derivative(
        (x_curr, y_curr) -> x_5(x_curr, y_curr; N=N, M=M),
        (x_curr, y_curr) -> dx_x_5(x_curr, y_curr; N=N, M=M),
        (x_curr, y_curr) -> dy_x_5(x_curr, y_curr; N=N, M=M),
        x,
        y;
        steps=steps,
    )
    println()

    println("x_6")
    test_xy_derivative(
        (x_curr, y_curr) -> x_6(x_curr, y_curr; N=N, M=M),
        (x_curr, y_curr) -> dx_x_6(x_curr, y_curr; N=N, M=M),
        (x_curr, y_curr) -> dy_x_6(x_curr, y_curr; N=N, M=M),
        x,
        y;
        steps=steps,
    )
end

"""
Вычисляет якобиан отображения `(s,x,y) -> ((1-s)x_4_prime + s x_4, x_5, x_6)`.
"""
function boundary_jacobian(x, y; N=N, M=M)
    return (x_4(x, y; N=N, M=M) - x_4_prime(x, y; N=N, M=M)) * (
        dx_x_5(x, y; N=N, M=M) * dy_x_6(x, y; N=N, M=M) -
        dy_x_5(x, y; N=N, M=M) * dx_x_6(x, y; N=N, M=M)
    )
end

"""
Вычисляет интегранд после аналитически дешёвого интегрирования по `s`.
"""
function boundary_integrand_xy(w, x, y; N=N, M=M, absolute_jacobian=false, atol=1e-9, rtol=1e-7)
    y45 = x_4_prime(x, y; N=N, M=M)
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

    jac = boundary_jacobian(x, y; N=N, M=M)
    if absolute_jacobian
        jac = abs(jac)
    end

    return inner_s * jac
end

"""
Вычисляет интеграл из последней формулы по заданным границам `x` и `y`.

Интегрирование выполняется в порядке `s -> y -> x`, чтобы внутренний интеграл по `s`
не обращался к тэта-функциям повторно.
"""
function integrate_boundary_volume(w, x_bounds, y_bounds; N=N, M=M, absolute_jacobian=false, atol=1e-9, rtol=1e-7)
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

    #test_der_x_on_surface_new(u, p, du, dp)
end
