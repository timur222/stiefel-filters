#??? Геометрия для сплошного изображения двух параметрических областей.
const _SolidPoint = NTuple{3, Float64}

_solid_add(a::_SolidPoint, b::_SolidPoint) = (a[1] + b[1], a[2] + b[2], a[3] + b[3]) #???
_solid_sub(a::_SolidPoint, b::_SolidPoint) = (a[1] - b[1], a[2] - b[2], a[3] - b[3]) #???
_solid_mul(c::Real, a::_SolidPoint) = (c * a[1], c * a[2], c * a[3]) #???
_solid_dot(a::_SolidPoint, b::_SolidPoint) = a[1] * b[1] + a[2] * b[2] + a[3] * b[3] #???
_solid_cross(a::_SolidPoint, b::_SolidPoint) = (a[2] * b[3] - a[3] * b[2],
                                                   a[3] * b[1] - a[1] * b[3],
                                                   a[1] * b[2] - a[2] * b[1]) #???
_solid_norm(a::_SolidPoint) = sqrt(_solid_dot(a, a)) #???
_solid_plane_value(p::_SolidPoint, normal::_SolidPoint, offset::Float64) =
    _solid_dot(p, normal) + offset #???

function _solid_near(a::_SolidPoint, b::_SolidPoint; tol=1e-8) #???
    return maximum(abs(a[k] - b[k]) for k in 1:3) <= tol
end

function _solid_clean_polygon(poly::Vector{_SolidPoint}) #???
    clean = _SolidPoint[]
    for p in poly
        (isempty(clean) || !_solid_near(p, clean[end])) && push!(clean, p)
    end
    length(clean) > 1 && _solid_near(clean[1], clean[end]) && pop!(clean)
    return clean
end

function _solid_clip_polygon(poly::Vector{_SolidPoint}, normal::_SolidPoint, offset::Float64) #???
    result = _SolidPoint[]
    for k in eachindex(poly)
        p, q = poly[k], poly[mod1(k + 1, length(poly))]
        vp = _solid_plane_value(p, normal, offset)
        vq = _solid_plane_value(q, normal, offset)
        inside_p, inside_q = vp >= -1e-10, vq >= -1e-10
        if inside_p != inside_q
            fraction = clamp(vp / (vp - vq), 0.0, 1.0)
            cut = _solid_add(p, _solid_mul(fraction, _solid_sub(q, p)))
            push!(result, cut)
        end
        inside_q && push!(result, q)
    end
    return _solid_clean_polygon(result)
end

function _solid_cap(points::Vector{_SolidPoint}, normal::_SolidPoint) #???
    unique_points = _SolidPoint[]
    for p in points
        all(q -> !_solid_near(p, q), unique_points) && push!(unique_points, p)
    end
    length(unique_points) >= 3 || return _SolidPoint[]

    center = _solid_mul(1 / length(unique_points), reduce(_solid_add, unique_points))
    reference = abs(normal[1]) < 0.9 ? (1.0, 0.0, 0.0) : (0.0, 1.0, 0.0)
    e1 = _solid_cross(normal, reference)
    e1 = _solid_mul(1 / _solid_norm(e1), e1)
    e2 = _solid_cross(normal, e1)
    sort!(unique_points; by=p -> begin
        delta = _solid_sub(p, center)
        atan(_solid_dot(delta, e2), _solid_dot(delta, e1))
    end)
    reverse!(unique_points)  # outward normal points opposite to the kept half-space
    return unique_points
end

#??? Пространство допустимых точек: X >= 0, Y >= X, Z >= Y, Z <= 100.
const _SOLID_PLANES = (
    ((1.0, 0.0, 0.0), 0.0),
    ((-1.0, 1.0, 0.0), 0.0),
    ((0.0, -1.0, 1.0), 0.0),
    ((0.0, 0.0, -1.0), 100.0),
)

function _solid_clip_tetra(tet::NTuple{4, _SolidPoint}) #???
    a, b, c, d = tet
    orientation = _solid_dot(_solid_cross(_solid_sub(b, a), _solid_sub(c, a)), _solid_sub(d, a))
    abs(orientation) > 1e-12 || return Vector{Vector{_SolidPoint}}()
    faces = Vector{Vector{_SolidPoint}}([
        [a, c, b], [a, b, d], [a, d, c], [b, c, d],
    ])
    orientation < 0 && foreach(reverse!, faces)

    for (normal, offset) in _SOLID_PLANES
        all_vertices = (p for face in faces for p in face)
        values = [_solid_plane_value(p, normal, offset) for p in all_vertices]
        crosses = any(v -> v > 1e-9, values) && any(v -> v < -1e-9, values)
        clipped_faces = Vector{Vector{_SolidPoint}}()
        cap_points = _SolidPoint[]
        for face in faces
            clipped = _solid_clip_polygon(face, normal, offset)
            if length(clipped) >= 3
                push!(clipped_faces, clipped)
                if crosses
                    append!(cap_points, (p for p in clipped if
                        abs(_solid_plane_value(p, normal, offset)) <= 1e-7))
                end
            end
        end
        if crosses
            cap = _solid_cap(cap_points, normal)
            length(cap) >= 3 && push!(clipped_faces, cap)
        end
        faces = clipped_faces
        isempty(faces) && break
    end
    return faces
end

function _solid_mesh_from_grid(lower::Matrix{_SolidPoint}, upper::Matrix{_SolidPoint}) #???
    nx, ny = size(lower)
    vertices = _SolidPoint[]
    visible_cells = 0
    #??? Freudenthal-разбиение каждой параметрической ячейки на шесть тетраэдров.
    tetrahedra = ((1, 2, 4, 8), (1, 2, 6, 8), (1, 5, 6, 8),
                  (1, 5, 7, 8), (1, 3, 7, 8), (1, 3, 4, 8))

    for j in 1:(ny - 1), i in 1:(nx - 1)
        cube = (
            lower[i, j], lower[i + 1, j], lower[i, j + 1], lower[i + 1, j + 1],
            upper[i, j], upper[i + 1, j], upper[i, j + 1], upper[i + 1, j + 1],
        )
        all(p -> all(isfinite, p), cube) || continue
        before = length(vertices)
        for ids in tetrahedra
            for face in _solid_clip_tetra(ntuple(k -> cube[ids[k]], 4))
                a = face[1]
                for k in 2:(length(face) - 1)
                    b, c = face[k], face[k + 1]
                    area2 = _solid_norm(_solid_cross(_solid_sub(b, a), _solid_sub(c, a)))
                    area2 > 1e-9 && append!(vertices, (a, b, c))
                end
            end
        end
        length(vertices) > before && (visible_cells += 1)
    end
    faces = permutedims(reshape(collect(1:length(vertices)), 3, length(vertices) ÷ 3))
    return vertices, faces, visible_cells
end

"""
    _draw_solid_volumes_for_M!(ax, region_color; N, M, x_bounds, y_bounds, nx, ny)

Добавляет к существующей оси две области для одного значения `M`.
"""
# !!! function plot_two_volumes(; N=8, M=4, x_bounds=(1e-3, 10.0),
# !!!                           y_bounds=(1e-3, 10.0), nx=31, ny=31)
function _draw_solid_volumes_for_M!(ax, region_color; N=8, M=4, #???
                                     x_bounds=(1e-3, 10.0),
                                     y_bounds=(1e-3, 10.0), nx=31, ny=31)
    0 < x_bounds[1] < x_bounds[2] || throw(ArgumentError("нужно 0 < x_bounds[1] < x_bounds[2]"))
    0 < y_bounds[1] < y_bounds[2] || throw(ArgumentError("нужно 0 < y_bounds[1] < y_bounds[2]"))
    nx >= 2 && ny >= 2 || throw(ArgumentError("nx и ny должны быть не меньше 2"))

    s_grid = exp.(range(log(x_bounds[1]), log(x_bounds[2]); length=nx))
    t_grid = exp.(range(log(y_bounds[1]), log(y_bounds[2]); length=ny))
    x4 = fill(NaN, nx, ny)
    y5 = fill(NaN, nx, ny)
    z6 = fill(NaN, nx, ny)
    for j in eachindex(t_grid), i in eachindex(s_grid)
        s, t = s_grid[i], t_grid[j]
        x4[i, j] = _plot_value(() -> x_4(s, t; N=N, M=M))
        y5[i, j] = _plot_value(() -> x_5(s, t; N=N, M=M))
        z6[i, j] = _plot_value(() -> x_6(s, t; N=N, M=M))
    end

    need_x4prime = falses(nx, ny)
    need_y5prime = falses(nx, ny)
    for j in 1:(ny - 1), i in 1:(nx - 1)
        corners = ((i, j), (i + 1, j), (i + 1, j + 1), (i, j + 1))
        all(((ii, jj),) -> isfinite(x4[ii, jj]) && isfinite(y5[ii, jj]) && isfinite(z6[ii, jj]),
            corners) || continue
        if maximum(y5[ii, jj] for (ii, jj) in corners) > 0 &&
           maximum(z6[ii, jj] - y5[ii, jj] for (ii, jj) in corners) > 0 &&
           minimum(z6[ii, jj] for (ii, jj) in corners) < 100
            foreach(((ii, jj),) -> (need_x4prime[ii, jj] = true), corners)
        end
        if maximum(x4[ii, jj] for (ii, jj) in corners) > 0 &&
           maximum(z6[ii, jj] - x4[ii, jj] for (ii, jj) in corners) > 0 &&
           minimum(z6[ii, jj] for (ii, jj) in corners) < 100
            foreach(((ii, jj),) -> (need_y5prime[ii, jj] = true), corners)
        end
    end

    x4prime = fill(NaN, nx, ny)
    y5prime = fill(NaN, nx, ny)
    for j in eachindex(t_grid), i in eachindex(s_grid)
        s, t = s_grid[i], t_grid[j]
        need_x4prime[i, j] &&
            (x4prime[i, j] = _plot_value(() -> x_4_prime(s, t; N=N, M=M)))
        need_y5prime[i, j] &&
            (y5prime[i, j] = _plot_value(() -> x_5_prime(s, t; N=N, M=M)))
    end

    invalid = (NaN, NaN, NaN)
    base = fill(invalid, nx, ny)
    end1 = fill(invalid, nx, ny)
    end2 = fill(invalid, nx, ny)
    for j in 1:ny, i in 1:nx
        if all(isfinite, (x4[i, j], y5[i, j], z6[i, j]))
            base[i, j] = (x4[i, j], y5[i, j], z6[i, j])
            isfinite(x4prime[i, j]) && (end1[i, j] = (x4prime[i, j], y5[i, j], z6[i, j]))
            isfinite(y5prime[i, j]) && (end2[i, j] = (x4[i, j], y5prime[i, j], z6[i, j]))
        end
    end

    raw1, faces1, cells1 = _solid_mesh_from_grid(base, end1)
    raw2, faces2, cells2 = _solid_mesh_from_grid(base, end2)
    # !!! @info "Непустые ячейки после обрезки" первая_область=cells1 вторая_область=cells2
    # !!! (cells1 == 0 || cells2 == 0) &&
    # !!!     @warn "Одна из областей не попала в сетку; измените диапазоны или увеличьте nx/ny."
    @info "Непустые ячейки после обрезки" M первая_область=cells1 вторая_область=cells2 #???

    # !!! fig = Makie.Figure(size=(1000, 800))
    # !!! ax = Makie.Axis3(fig[1, 1]; xlabel="X", ylabel="Y", zlabel="Z",
    # !!!                  title="Область между x₄, x₄′ и x₅, x₅′")
    # !!! region_color = Makie.RGBAf(0.12, 0.56, 1.0, 1.0)
    cells1 > 0 && Makie.mesh!(ax, [Makie.Point3f(p...) for p in raw1], faces1;
                              color=region_color, shading=Makie.NoShading)
    cells2 > 0 && Makie.mesh!(ax, [Makie.Point3f(p...) for p in raw2], faces2;
                              color=region_color, shading=Makie.NoShading)
    # !!! Makie.limits!(ax, 0, 100, 0, 100, 0, 100)
    # !!! display(fig)
    # !!! return fig
    return cells1, cells2 #???
end

"""
    plot_two_volumes(; N=8, M=1:(N-1), x_bounds=(1e-3, 10.0),
                       y_bounds=(1e-3, 10.0), nx=31, ny=31)

Рисует на одном графике обе сплошные области для каждого `M` из `1:N-1`.
Можно передать одно целое `M`, диапазон или список значений `M`.
"""
function plot_two_volumes(; N=8, M=1:(N-1), x_bounds=(1e-3, 10.0),
                          y_bounds=(1e-3, 10.0), nx=31, ny=31) #???
    N isa Integer && N >= 2 || throw(ArgumentError("N должно быть целым числом не меньше 2"))
    M_values = M isa Integer ? [M] : collect(M)
    isempty(M_values) && throw(ArgumentError("список M не может быть пустым"))
    all(m -> m isa Integer && 1 <= m < N, M_values) ||
        throw(ArgumentError("каждое M должно быть целым числом из 1:N-1"))
    M_values = unique(Int.(M_values))

    fig = Makie.Figure(size=(1000, 800))
    ax = Makie.Axis3(fig[1, 1]; xlabel="X", ylabel="Y", zlabel="Z",
                     title="Области для M = $(join(M_values, ", "))")
    region_color = Makie.RGBAf(0.12, 0.56, 1.0, 1.0)
    total_cells = 0
    for m in M_values
        cells1, cells2 = _draw_solid_volumes_for_M!(ax, region_color; N=N, M=m,
            x_bounds=x_bounds, y_bounds=y_bounds, nx=nx, ny=ny)
        total_cells += cells1 + cells2
    end
    total_cells == 0 && @warn "Ни одна область не попала в сетку; измените диапазоны или увеличьте nx/ny."
    Makie.limits!(ax, 0, 100, 0, 100, 0, 100)
    display(fig)
    return fig
end
