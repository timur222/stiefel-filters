#??? Две параметрические поверхности внутри области 0 < X < Y < Z < 100.
include(joinpath(@__DIR__, "../stiefel/points.jl"))

function _plot_value(f) #???
    try
        value = f()
        return isfinite(value) ? Float64(value) : NaN
    catch err
        if err isa DomainError || err isa DivideError || err isa OverflowError ||
           err isa LinearAlgebra.SingularException || err isa LinearAlgebra.PosDefException
            return NaN
        end
        rethrow()
    end
end

#??? Пересечение отрезка a + l(b-a), 0 <= l <= 1, с интервалом low < value < high.
function _clipped_l_range(a, b, low, high) #???
    all(isfinite, (a, b, low, high)) && low < high || return nothing
    if a == b
        return low < a < high ? (0.0, 1.0) : nothing
    end
    l1, l2 = (low - a) / (b - a), (high - a) / (b - a)
    l_min = max(0.0, min(l1, l2))
    l_max = min(1.0, max(l1, l2))
    return l_min < l_max ? (l_min, l_max) : nothing
end

function _volume_boundary_mesh(lower, upper) #???
    nx, ny = size(lower)
    valid = falses(nx - 1, ny - 1)
    for j in 1:(ny - 1), i in 1:(nx - 1)
        valid[i, j] = all(
            all(isfinite, lower[ii, jj]) && all(isfinite, upper[ii, jj])
            for (ii, jj) in ((i, j), (i + 1, j), (i + 1, j + 1), (i, j + 1))
        )
    end

    vertices = Makie.Point3f[]
    indices = Int[]
    function add_quad(a, b, c, d)
        offset = length(vertices)
        append!(vertices, (Makie.Point3f(a...), Makie.Point3f(b...),
                           Makie.Point3f(c...), Makie.Point3f(d...)))
        append!(indices, (offset + 1, offset + 2, offset + 3,
                          offset + 1, offset + 3, offset + 4))
    end

    #??? Рисуем только внешние грани ячеек: внутренние грани не затемняют объём.
    for j in 1:(ny - 1), i in 1:(nx - 1)
        valid[i, j] || continue
        a, b, c, d = lower[i, j], lower[i + 1, j],
                     lower[i + 1, j + 1], lower[i, j + 1]
        e, f, g, h = upper[i, j], upper[i + 1, j],
                     upper[i + 1, j + 1], upper[i, j + 1]
        add_quad(a, d, c, b)  # l = l_min
        add_quad(e, f, g, h)  # l = l_max
        (i == 1 || !valid[i - 1, j]) && add_quad(a, e, h, d)
        (i == nx - 1 || !valid[i + 1, j]) && add_quad(b, c, g, f)
        (j == 1 || !valid[i, j - 1]) && add_quad(a, b, f, e)
        (j == ny - 1 || !valid[i, j + 1]) && add_quad(d, h, g, c)
    end

    faces = permutedims(reshape(indices, 3, length(indices) ÷ 3))
    return vertices, faces, count(valid)
end

#= !!! Прежний способ построения объёмов оставлен для сравнения. Он выбрасывал граничные ячейки.
"""
    plot_two_volumes(; N=8, M=4, x_bounds=(1e-3, 10.0),
                       y_bounds=(1e-3, 10.0), nx=31, ny=31)

Рисует две области, полученные при 0 < l < 1 из
`(l*x_4_prime(s,t) + (1-l)*x_4(s,t), x_5(s,t), x_6(s,t))` и
`(x_4(s,t), l*x_5_prime(s,t) + (1-l)*x_5(s,t), x_6(s,t))`.
Показывает только их части внутри 0 < X < Y < Z < 100.
"""
function plot_two_volumes(; N=8, M=4, x_bounds=(1e-3, 10.0),
                          y_bounds=(1e-3, 10.0), nx=31, ny=31) #???
    0 < x_bounds[1] < x_bounds[2] || throw(ArgumentError("нужно 0 < x_bounds[1] < x_bounds[2]"))
    0 < y_bounds[1] < y_bounds[2] || throw(ArgumentError("нужно 0 < y_bounds[1] < y_bounds[2]"))
    nx >= 2 && ny >= 2 || throw(ArgumentError("nx и ny должны быть не меньше 2"))

    s_grid = exp.(range(log(x_bounds[1]), log(x_bounds[2]); length=nx))
    t_grid = exp.(range(log(y_bounds[1]), log(y_bounds[2]); length=ny))
    missing_point = (NaN, NaN, NaN)
    lower1, upper1 = fill(missing_point, nx, ny), fill(missing_point, nx, ny)
    lower2, upper2 = fill(missing_point, nx, ny), fill(missing_point, nx, ny)

    for j in eachindex(t_grid), i in eachindex(s_grid)
        s, t = s_grid[i], t_grid[j]
        z6 = _plot_value(() -> x_6(s, t; N=N, M=M))
        0 < z6 < 100 || continue

        x4 = _plot_value(() -> x_4(s, t; N=N, M=M))
        y5 = _plot_value(() -> x_5(s, t; N=N, M=M))

        if 0 < y5 < z6
            x4prime = _plot_value(() -> x_4_prime(s, t; N=N, M=M))
            interval = _clipped_l_range(x4, x4prime, 0.0, y5)
            if interval !== nothing
                l_min, l_max = interval
                lower1[i, j] = (x4 + l_min * (x4prime - x4), y5, z6)
                upper1[i, j] = (x4 + l_max * (x4prime - x4), y5, z6)
            end
        end

        if 0 < x4 < z6
            y5prime = _plot_value(() -> x_5_prime(s, t; N=N, M=M))
            interval = _clipped_l_range(y5, y5prime, x4, z6)
            if interval !== nothing
                l_min, l_max = interval
                lower2[i, j] = (x4, y5 + l_min * (y5prime - y5), z6)
                upper2[i, j] = (x4, y5 + l_max * (y5prime - y5), z6)
            end
        end
    end

    vertices1, faces1, cells1 = _volume_boundary_mesh(lower1, upper1)
    vertices2, faces2, cells2 = _volume_boundary_mesh(lower2, upper2)
    @info "Ячейки внутри 0 < X < Y < Z < 100" первая_область=cells1 вторая_область=cells2
    if cells1 == 0 || cells2 == 0
        @warn "Одна из областей не содержит целых ячеек сетки; увеличьте nx/ny или измените диапазоны s,t."
    end

    #??? Прозрачные замкнутые оболочки показывают объём каждой области.
    fig = Makie.Figure(size=(1000, 800))
    # !!! ax = Makie.Axis3(fig[1, 1]; xlabel="X", ylabel="Y", zlabel="Z",
    # !!!                  title="Синяя: интервал x₄→x₄′; оранжевая: интервал x₅→x₅′")
    ax = Makie.Axis3(fig[1, 1]; xlabel="X", ylabel="Y", zlabel="Z", #???
                     title="Область между x₄, x₄′ и x₅, x₅′") #???
    region_color = Makie.RGBAf(0.12, 0.56, 1.0, 0.35) #???
    if cells1 > 0
        # !!! Makie.mesh!(ax, vertices1, faces1; color=Makie.RGBAf(0.12, 0.56, 1.0, 0.35),
        # !!!             transparency=true, shading=Makie.NoShading)
        Makie.mesh!(ax, vertices1, faces1; color=region_color, #???
                    transparency=true, shading=Makie.NoShading) #???
    end
    if cells2 > 0
        # !!! Makie.mesh!(ax, vertices2, faces2; color=Makie.RGBAf(1.0, 0.65, 0.0, 0.35),
        # !!!             transparency=true, shading=Makie.NoShading)
        Makie.mesh!(ax, vertices2, faces2; color=region_color, #???
                    transparency=true, shading=Makie.NoShading) #???
    end
    Makie.limits!(ax, 0, 100, 0, 100, 0, 100)
    display(fig)
    return fig
end
=#

"""
    plot_two_surfaces(; N=8, M=4, x_bounds=(1e-3, 10.0),
                        y_bounds=(1e-3, 10.0), nx=31, ny=31)

Рисует на одном графике поверхности
`(x_4_prime(u,v), x_5(u,v), x_6(u,v))` и
`(x_4(u,v), x_5_prime(u,v), x_6(u,v))` при `u,v > 0`.
Оставляет только точки с `0 < X < Y < Z < 100`.
`x_bounds` и `y_bounds` задают диапазоны входных параметров, а не координат X,Y,Z.
"""
function plot_two_surfaces(; N=8, M=4, x_bounds=(1e-3, 10.0),
                           y_bounds=(1e-3, 10.0), nx=31, ny=31) #???
    0 < x_bounds[1] < x_bounds[2] || throw(ArgumentError("нужно 0 < x_bounds[1] < x_bounds[2]"))
    0 < y_bounds[1] < y_bounds[2] || throw(ArgumentError("нужно 0 < y_bounds[1] < y_bounds[2]"))
    nx >= 2 && ny >= 2 || throw(ArgumentError("nx и ny должны быть не меньше 2"))

    #??? Логарифмическая сетка захватывает малые и большие положительные параметры.
    u_grid = exp.(range(log(x_bounds[1]), log(x_bounds[2]); length=nx))
    v_grid = exp.(range(log(y_bounds[1]), log(y_bounds[2]); length=ny))
    X1, Y1, Z1 = (fill(NaN, nx, ny) for _ in 1:3)
    X2, Y2, Z2 = (fill(NaN, nx, ny) for _ in 1:3)

    for j in eachindex(v_grid), i in eachindex(u_grid)
        u, v = u_grid[i], v_grid[j]
        z6 = _plot_value(() -> x_6(u, v; N=N, M=M))
        0 < z6 < 100 || continue

        y5 = _plot_value(() -> x_5(u, v; N=N, M=M))
        if 0 < y5 < z6
            x4prime = _plot_value(() -> x_4_prime(u, v; N=N, M=M))
            if 0 < x4prime < y5
                X1[i, j], Y1[i, j], Z1[i, j] = x4prime, y5, z6
            end
        end

        x4 = _plot_value(() -> x_4(u, v; N=N, M=M))
        if 0 < x4 < z6
            y5prime = _plot_value(() -> x_5_prime(u, v; N=N, M=M))
            if x4 < y5prime < z6
                X2[i, j], Y2[i, j], Z2[i, j] = x4, y5prime, z6
            end
        end
    end

    n1, n2 = count(isfinite, X1), count(isfinite, X2)
    @info "Точки в области 0 < X < Y < Z < 100" первая_поверхность=n1 вторая_поверхность=n2
    if n1 == 0 && n2 == 0
        @warn "Обе поверхности пусты в выбранном диапазоне входных параметров; измените x_bounds/y_bounds."
    end

    #??? NaN удаляет прилегающие треугольники за пределами заданной области.
    fig = Makie.Figure(size=(1000, 800))
    ax = Makie.Axis3(fig[1, 1]; xlabel="X", ylabel="Y", zlabel="Z",
                     title="Синяя: (x₄′, x₅, x₆); оранжевая: (x₄, x₅′, x₆)")
    # !!! Makie.surface!(ax, X1, Y1, Z1; color=(:dodgerblue, 0.65), transparency=true)
    # !!! Makie.surface!(ax, X2, Y2, Z2; color=(:orange, 0.65), transparency=true)
    color1 = fill(Makie.RGBAf(0.12, 0.56, 1.0, 0.65), size(X1)) #???
    color2 = fill(Makie.RGBAf(1.0, 0.65, 0.0, 0.65), size(X2)) #???
    Makie.surface!(ax, X1, Y1, Z1; color=color1, transparency=true) #???
    Makie.surface!(ax, X2, Y2, Z2; color=color2, transparency=true) #???
    Makie.limits!(ax, 0, 100, 0, 100, 0, 100)
    display(fig)
    return fig
end

include(joinpath(@__DIR__, "solid_volume_mesh.jl")) #???
include(joinpath(@__DIR__, "volume_slice.jl")) #???
