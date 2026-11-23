#??? Сечение объединения параметрических областей плоскостью Z = z0.
const _SlicePoint = NTuple{2, Float64}

function _slice_crosses(z, corners, z0) #???
    values = (z[i, j] for (i, j) in corners)
    all(isfinite, values) || return false
    levels = [z[i, j] for (i, j) in corners]
    return minimum(levels) <= z0 <= maximum(levels) && minimum(levels) < maximum(levels)
end

function _slice_edge_value(a, b, za, zb, z0) #???
    if (za <= z0 < zb) || (zb <= z0 < za)
        weight = (z0 - za) / (zb - za)
        return ntuple(k -> a[k] + weight * (b[k] - a[k]), 4)
    end
    return nothing
end

# !!! function _slice_segment(corners, z, data, z0)
function _slice_segment(corners, z, data, z0; refine_edge=nothing) #???
    cut = NTuple{4, Float64}[]
    for (k1, k2) in ((1, 2), (2, 3), (3, 1))
        i1, j1 = corners[k1]
        i2, j2 = corners[k2]
        # !!! point = _slice_edge_value(data[i1, j1], data[i2, j2], z[i1, j1], z[i2, j2], z0)
        point = isnothing(refine_edge) ? #???
            _slice_edge_value(data[i1, j1], data[i2, j2], z[i1, j1], z[i2, j2], z0) : #???
            refine_edge(i1, j1, i2, j2) #???
        point !== nothing && push!(cut, point)
    end
    return length(cut) == 2 ? (cut[1], cut[2]) : nothing
end

# !!! Уточняли пересечение методом секущих без производных x₆.
#??? Уточняем пересечение уровнем x₆=z0 методом Ньютона внутри скобки.
#??? При выходе шага за ребро или вырождении производной делим скобку пополам.
function _slice_refined_edge_value(i1, j1, i2, j2, za, zb, z0, s_grid, t_grid;
                                   # !!! N, M, maxiter=30, rtol=1e-9) #???
                                   N, M, maxiter=10, rtol=1e-9) #??? Ограничиваем цену безуспешной попытки Ньютона.
    (za <= z0 < zb || zb <= z0 < za) || return nothing #???
    s1, s2 = s_grid[i1], s_grid[i2] #???
    t1, t2 = t_grid[j1], t_grid[j2] #???
    ds, dt = s2 - s1, t2 - t1 #???
    low, high = 0.0, 1.0 #???
    # !!! flow, fhigh = za - z0, zb - z0 #??? Для метода секущих.
    flow = za - z0 #???
    λ = clamp((z0 - za) / (zb - za), low, high) #???
    tolerance = rtol * max(1.0, abs(z0)) #???
    # !!! last_side = 0 #??? Ускорение Illinois для метода секущих.
    for _ in 1:maxiter #???
        s, t = s1 + λ * ds, t1 + λ * dt #???
        height = _plot_value(() -> x_6(s, t; N=N, M=M)) #???
        isfinite(height) || return nothing #???
        residual = height - z0 #???
        if abs(residual) <= tolerance #???
            values = ntuple(k -> _plot_value(() -> #???
                (x_4, x_5, x_4_prime, x_5_prime)[k](s, t; N=N, M=M)), 4) #???
            return all(isfinite, values) ? values : nothing #???
        end #???
        if (residual < 0) == (flow < 0) #???
            low, flow = λ, residual #???
            # !!! last_side == 1 && (fhigh *= 0.5) #???
            # !!! last_side = 1 #???
        else #???
            # !!! high, fhigh = λ, residual #???
            high = λ #???
            # !!! last_side == -1 && (flow *= 0.5) #???
            # !!! last_side = -1 #???
        end #???
        derivative = _plot_value(() -> dx_x_6(s, t; N=N, M=M)) * ds + #???
                     _plot_value(() -> dy_x_6(s, t; N=N, M=M)) * dt #???
        # !!! span = high - low #???
        # !!! weight = abs(flow) / (abs(flow) + abs(fhigh)) #???
        # !!! candidate = low + span * weight #???
        candidate = isfinite(derivative) && derivative != 0 ? λ - residual / derivative : NaN #???
        λ = isfinite(candidate) && low < candidate < high ? candidate : (low + high) / 2 #???
    end #???
    return nothing #???
end #???

function _slice_clip_polygon(poly::Vector{_SlicePoint}, value) #???
    result = _SlicePoint[]
    for k in eachindex(poly)
        p, q = poly[k], poly[mod1(k + 1, length(poly))]
        vp, vq = value(p), value(q)
        inside_p, inside_q = vp >= -1e-10, vq >= -1e-10
        if inside_p != inside_q
            weight = clamp(vp / (vp - vq), 0.0, 1.0)
            push!(result, (p[1] + weight * (q[1] - p[1]),
                           p[2] + weight * (q[2] - p[2])))
        end
        inside_q && push!(result, q)
    end
    return result
end

function _slice_domain_clip(poly::Vector{_SlicePoint}, z0) #???
    all(p -> all(isfinite, p), poly) || return _SlicePoint[]
    for value in (p -> p[1], p -> p[2] - p[1], p -> z0 - p[2])
        poly = _slice_clip_polygon(poly, value)
        length(poly) >= 3 || return _SlicePoint[]
    end
    return poly
end

function _slice_area(poly::Vector{_SlicePoint}) #???
    total = 0.0
    for k in eachindex(poly)
        p, q = poly[k], poly[mod1(k + 1, length(poly))]
        total += p[1] * q[2] - q[1] * p[2]
    end
    return abs(total) / 2
end

function _slice_quads(a, b, kind) #???
    all(isfinite, (a[1], a[2], b[1], b[2])) || return Vector{Vector{_SlicePoint}}()
    endpoint = kind == 1 ? 3 : 4
    all(isfinite, (a[endpoint], b[endpoint])) || return Vector{Vector{_SlicePoint}}()
    displacement_a = a[endpoint] - a[kind]
    displacement_b = b[endpoint] - b[kind]
    segments = ((a, b),)
    if displacement_a * displacement_b < 0
        weight = displacement_a / (displacement_a - displacement_b)
        middle = ntuple(k -> a[k] + weight * (b[k] - a[k]), 4)
        segments = ((a, middle), (middle, b))
    end
    polygons = Vector{Vector{_SlicePoint}}()
    for (left, right) in segments
        if kind == 1
            poly = _SlicePoint[(left[1], left[2]), (right[1], right[2]),
                               (right[3], right[2]), (left[3], left[2])]
        else
            poly = _SlicePoint[(left[1], left[2]), (right[1], right[2]),
                               (right[1], right[4]), (left[1], left[4])]
        end
        push!(polygons, poly)
    end
    return polygons
end

function _slice_contains(poly::Vector{_SlicePoint}, x, y) #???
    inside = false
    for k in eachindex(poly)
        a, b = poly[k], poly[mod1(k + 1, length(poly))]
        if (a[2] > y) != (b[2] > y)
            crossing_x = a[1] + (y - a[2]) * (b[1] - a[1]) / (b[2] - a[2])
            x < crossing_x && (inside = !inside)
        end
    end
    return inside
end

function _slice_paint!(mask, coordinates, poly) #???
    _slice_area(poly) > 1e-10 || return false
    x_min, x_max = extrema(p[1] for p in poly)
    y_min, y_max = extrema(p[2] for p in poly)
    i_min = searchsortedfirst(coordinates, x_min)
    i_max = searchsortedlast(coordinates, x_max)
    j_min = searchsortedfirst(coordinates, y_min)
    j_max = searchsortedlast(coordinates, y_max)
    for j in j_min:j_max, i in i_min:i_max
        mask[i, j] || (mask[i, j] = _slice_contains(poly, coordinates[i], coordinates[j]))
    end
    return true
end

#= !!! Возвращали вместе с координатой найденный корень для продолжения ветви на соседнем узле.
function _slice_projected_value(f) #???
    root = Ref{Union{Nothing, ComplexF64}}(nothing) #???
    value = _plot_value(() -> begin #???
        projected, next_root = f() #???
        root[] = next_root #???
        projected #???
    end) #???
    return value, isfinite(value) ? root[] : nothing #???
end #???
=#


# !!! function _slice_mask_for_M(z0; N, M, x_bounds, y_bounds, nx, ny, pixels) #???
# !!! function _slice_mask_for_M(z0; N, M, x_bounds, y_bounds, nx, ny, pixels, coordinates=nothing)
function _slice_mask_for_M(z0; N, M, x_bounds, y_bounds, nx, ny, pixels,
                           coordinates=nothing, refine_intersections=true) #???
    stage_started = time_ns() #??? Временные отметки покажут, на каком этапе задерживается расчёт.
    progress_step = max(1, cld(ny, 5)) #???
    @info "Сечение: вычисление x₆" M=M узлов=nx * ny #???
    s_grid = exp.(range(log(x_bounds[1]), log(x_bounds[2]); length=nx))
    t_grid = exp.(range(log(y_bounds[1]), log(y_bounds[2]); length=ny))
    z6 = fill(NaN, nx, ny)
    for j in 1:ny, i in 1:nx
        z6[i, j] = _plot_value(() -> x_6(s_grid[i], t_grid[j]; N=N, M=M))
        if i == nx && (j % progress_step == 0 || j == ny) #???
            @info "Сечение: x₆" M=M строк=j всего=ny секунд=round((time_ns() - stage_started) / 1e9; digits=1) #???
        end #???
    end

    needed = falses(nx, ny)
    for j in 1:(ny - 1), i in 1:(nx - 1)
        for corners in (((i, j), (i + 1, j), (i + 1, j + 1)),
                        ((i, j), (i + 1, j + 1), (i, j + 1)))
            _slice_crosses(z6, corners, z0) || continue
            for (ii, jj) in corners
                needed[ii, jj] = true
            end
        end
    end

    invalid = (NaN, NaN, NaN, NaN)
    data = fill(invalid, nx, ny)
    stage_started = time_ns() #???
    @info "Сечение: вычисление границ" M=M нужных_узлов=count(needed) #???
    #??? Возвращён расчёт каждой проекции независимо, как до предыдущей итерации.
    for j in 1:ny #???
        for i in 1:nx #???
            needed[i, j] || continue
            s, t = s_grid[i], t_grid[j]
            data[i, j] = (
                _plot_value(() -> x_4(s, t; N=N, M=M)),
                _plot_value(() -> x_5(s, t; N=N, M=M)),
                _plot_value(() -> x_4_prime(s, t; N=N, M=M)),
                _plot_value(() -> x_5_prime(s, t; N=N, M=M)),
            )
        end #???
        if j % progress_step == 0 || j == ny #???
            @info "Сечение: границы" M=M строк=j всего=ny секунд=round((time_ns() - stage_started) / 1e9; digits=1) #???
        end #???
    end
    @info "Сечение: границы вычислены" M=M секунд=round((time_ns() - stage_started) / 1e9; digits=1) #???
    #= !!! Продолжение корня и диагностика скачков из предыдущей итерации.
    x4_seeds = Matrix{Union{Nothing, ComplexF64}}(undef, nx, ny) #???
    x5_seeds = Matrix{Union{Nothing, ComplexF64}}(undef, nx, ny) #???
    fill!(x4_seeds, nothing) #???
    fill!(x5_seeds, nothing) #???
    fallback_projections = 0 #???
    max_root_jump = 0.0 #???
    for j in 1:ny, i in 1:nx #???
        needed[i, j] || continue #???
        s, t = s_grid[i], t_grid[j] #???
        seed4 = i > 1 ? x4_seeds[i - 1, j] : nothing #???
        seed5 = i > 1 ? x5_seeds[i - 1, j] : nothing #???
        seed4 === nothing && j > 1 && (seed4 = x4_seeds[i, j - 1]) #???
        seed5 === nothing && j > 1 && (seed5 = x5_seeds[i, j - 1]) #???
        x4p, next4 = _slice_projected_value(() -> x_4_prime( #???
            s, t; N=N, M=M, active_seed=seed4, return_seed=true)) #???
        x5p, next5 = _slice_projected_value(() -> x_5_prime( #???
            s, t; N=N, M=M, active_seed=seed5, return_seed=true)) #???
        if !isfinite(x4p) && seed4 !== nothing #???
            x4p, next4 = _slice_projected_value(() -> x_4_prime( #???
                s, t; N=N, M=M, return_seed=true)) #???
            fallback_projections += 1 #???
        end #???
        if !isfinite(x5p) && seed5 !== nothing #???
            x5p, next5 = _slice_projected_value(() -> x_5_prime( #???
                s, t; N=N, M=M, return_seed=true)) #???
            fallback_projections += 1 #???
        end #???
        x4_seeds[i, j], x5_seeds[i, j] = next4, next5 #???
        next4 !== nothing && seed4 !== nothing && #???
            (max_root_jump = max(max_root_jump, abs(next4 - seed4))) #???
        next5 !== nothing && seed5 !== nothing && #???
            (max_root_jump = max(max_root_jump, abs(next5 - seed5))) #???
        data[i, j] = ( #???
            _plot_value(() -> x_4(s, t; N=N, M=M)), #???
            _plot_value(() -> x_5(s, t; N=N, M=M)), #???
            x4p, x5p, #???
        ) #???
    end #???
    @info "Продолжение корня на сечении" M=M повторных_попыток=fallback_projections максимальный_скачок_корня=max_root_jump #???
    #??? Ищем самый большой скачок между соседними значениями границы до интерполяции многоугольников.
    max_boundary_jump = 0.0 #???
    jump_details = nothing #???
    boundary_names = (:x4, :x5, :x4_prime, :x5_prime) #???
    for j in 1:ny, i in 1:nx #???
        needed[i, j] || continue #???
        for (ii, jj) in ((i - 1, j), (i, j - 1)) #???
            ii >= 1 && jj >= 1 && needed[ii, jj] || continue #???
            for k in 1:4 #???
                current, previous = data[i, j][k], data[ii, jj][k] #???
                isfinite(current) && isfinite(previous) || continue #???
                jump = abs(current - previous) #???
                if jump > max_boundary_jump #???
                    max_boundary_jump = jump #???
                    jump_details = (граница=boundary_names[k], s=s_grid[i], t=t_grid[j], #???
                                    предыдущее=previous, текущее=current) #???
                end #???
            end #???
        end #???
    end #???
    @info "Наибольший скачок исходной границы" M=M величина=max_boundary_jump данные=jump_details #???
    =#
    #??? Считаем отброшенные проекции: они не должны порождать граничные многоугольники.
    # !!! invalid_projections = count(needed[i, j] && (!isfinite(data[i, j][3]) || !isfinite(data[i, j][4])) for j in 1:ny, i in 1:nx)
    invalid_x4 = count(needed[i, j] && !isfinite(data[i, j][3]) for j in 1:ny, i in 1:nx) #???
    invalid_x5 = count(needed[i, j] && !isfinite(data[i, j][4]) for j in 1:ny, i in 1:nx) #???
    invalid_projections = count(needed[i, j] && #???
        (!isfinite(data[i, j][3]) || !isfinite(data[i, j][4])) #???
        for j in 1:ny, i in 1:nx) #???
    invalid_projections > 0 && @info "Отброшенные проекции с невязкой theta" M=M точек=invalid_projections x₄prime=invalid_x4 x₅prime=invalid_x5 всего_нужных=count(needed) #???

    # !!! coordinates = range(0.0, Float64(z0); length=pixels)
    coordinates = isnothing(coordinates) ? range(0.0, Float64(z0); length=pixels) : coordinates #???
    mask = falses(pixels, pixels)
    polygons = 0
    edge_cache = Dict{Tuple{Tuple{Int, Int}, Tuple{Int, Int}}, Union{Nothing, NTuple{4, Float64}}}() #???
    refined_edges = Ref(0) #???
    fallback_edges = Ref(0) #??? Запасная линейная оценка сохраняет граничные многоугольники.
    rejected_edges = Ref(0) #???
    skip_refinement = Ref(false) #??? После серии неудач не повторяем дорогой поиск для каждого ребра.
    consecutive_failures = Ref(0) #??? Считаем неудачи подряд даже после нескольких успешных рёбер.
    stage_started = time_ns() #???
    @info "Сечение: уточнение пересечений и заливка" M=M #???
    function refine_edge(i1, j1, i2, j2) #???
        key = (i1, j1) < (i2, j2) ? ((i1, j1), (i2, j2)) : ((i2, j2), (i1, j1)) #???
        return get!(edge_cache, key) do #???
            za, zb = z6[i1, j1], z6[i2, j2] #???
            a, b = data[i1, j1], data[i2, j2] #???
            linear = _slice_edge_value(a, b, za, zb, z0) #???
            isnothing(linear) && return nothing #???
            all(isfinite, linear) || return nothing #??? Невалидные концы не могут дать корректный многоугольник.
            if skip_refinement[] #???
                fallback_edges[] += 1 #???
                return linear #???
            end #???
            refined = _slice_refined_edge_value(i1, j1, i2, j2, za, zb, z0,
                                                 s_grid, t_grid; N=N, M=M) #???
            if isnothing(refined) #???
                # !!! return nothing #??? Прежде отбрасывали всё ребро и теряли тёмную область.
                rejected_edges[] += 1 #???
                consecutive_failures[] += 1 #???
                fallback_edges[] += 1 #???
                # !!! if rejected_edges[] >= 64 && refined_edges[] == 0 #???
                if consecutive_failures[] >= 8 #???
                    skip_refinement[] = true #???
                    @warn "Уточнение x₆=z₀ не сходится; остальные рёбра оцениваются линейно" M=M подряд_неудач=consecutive_failures[] всего_проверено=refined_edges[] + rejected_edges[] #???
                end #???
                return linear #???
            end #???
            refined_edges[] += 1 #???
            consecutive_failures[] = 0 #???
            return refined #???
        end #???
    end #???
    # !!! corner_polygons = 0 #???
    # !!! boundary_tol = 1e-8 * max(1.0, abs(z0)) #???
    for j in 1:(ny - 1), i in 1:(nx - 1)
        for corners in (((i, j), (i + 1, j), (i + 1, j + 1)),
                        ((i, j), (i + 1, j + 1), (i, j + 1)))
            _slice_crosses(z6, corners, z0) || continue
            # !!! segment = _slice_segment(corners, z6, data, z0)
            segment = _slice_segment(corners, z6, data, z0; #???
                                     refine_edge=refine_intersections ? refine_edge : nothing) #???
            segment === nothing && continue
            for kind in 1:2, quad in _slice_quads(segment[1], segment[2], kind)
                poly = _slice_domain_clip(quad, z0)
                isempty(poly) && continue
                #= !!! Прежний фильтр удалял целые многоугольники, касавшиеся сразу X=0 и Y=Z.
                if minimum(point[1] for point in poly) <= boundary_tol && #???
                   maximum(point[2] for point in poly) >= z0 - boundary_tol #???
                    corner_polygons += 1 #???
                    continue #???
                end #???
                =#
                polygons += _slice_paint!(mask, coordinates, poly)
            end
        end
        if i == nx - 1 && (j % progress_step == 0 || j == ny - 1) #???
            @info "Сечение: пересечения" M=M строк=j всего=ny - 1 уточнено=refined_edges[] линейно=fallback_edges[] отклонено=rejected_edges[] секунд=round((time_ns() - stage_started) / 1e9; digits=1) #???
        end #???
    end
    # !!! corner_polygons > 0 && @info "Отброшенные многоугольники в углу X=0, Y=Z" M=M многоугольников=corner_polygons #???
    refine_intersections && @info "Уточнение пересечений x₆=z0" M=M уточнено=refined_edges[] линейно=fallback_edges[] отклонено=rejected_edges[] #???
    @info "Сечение: маска готова" M=M секунд=round((time_ns() - stage_started) / 1e9; digits=1) #???
    return mask, polygons
end

#??? Находим τ, при котором верхняя правая вершина прямоугольника имеет высоту Z=z0.
function _slice_gen1_tau(z0, N; rtol=1e-7, max_steps=60) #???
    height(τ) = x_gen_1(τ, N; special_point=4) #???
    lower = 1.0 / N #???
    upper = lower #???
    lower_height = height(lower) #???
    upper_height = lower_height #???
    if lower_height == z0 #???
        return lower #???
    elseif lower_height < z0 #???
        for _ in 1:max_steps #???
            lower, lower_height = upper, upper_height #???
            upper *= 2 #???
            upper_height = height(upper) #???
            upper_height >= z0 && break #???
        end #???
    else #???
        for _ in 1:max_steps #???
            upper, upper_height = lower, lower_height #???
            lower /= 2 #???
            lower_height = height(lower) #???
            lower_height <= z0 && break #???
        end #???
    end #???
    if !(lower_height <= z0 && z0 <= upper_height) #???
        throw(ArgumentError("Не удалось найти τ для сечения Z=$z0")) #???
    end #???
    for _ in 1:max_steps #???
        middle = (lower + upper) / 2 #???
        middle_height = height(middle) #???
        abs(middle_height - z0) <= rtol * max(1.0, z0) && return middle #???
        if middle_height < z0 #???
            lower = middle #???
        else #???
            upper = middle #???
        end #???
    end #???
    return (lower + upper) / 2 #???
end #???

#??? Два треугольника для каждого M в плоскости Z=z0, как в параметризации через s,t.
# !!! function _slice_gen1_masks(z0; N, M_values, pixels, coordinates=nothing, keep_masks=true)
function _slice_gen1_masks(z0; N, M_values, pixels, coordinates=nothing) #???
    τ = _slice_gen1_tau(z0, N) #???
    right = Float64[x_gen_1(τ, N; right_point=k) for k in 0:N] #???
    # !!! coordinates = range(0.0, Float64(z0); length=pixels) #???
    coordinates = isnothing(coordinates) ? range(0.0, Float64(z0); length=pixels) : coordinates #???
    masks = BitMatrix[] #???
    combined = falses(pixels, pixels) #???
    for m in M_values #???
        mask = falses(pixels, pixels) #???
        polygons = 0 #???
        for (a, b) in ((right[m], right[m + 1]), (right[m + 1], right[m + 2])) #???
            triangle = _SlicePoint[(a, a), (a, b), (b, b)] #???
            clipped = _slice_domain_clip(triangle, z0) #???
            isempty(clipped) && continue #???
            polygons += _slice_paint!(mask, coordinates, clipped) #???
        end #???
        @info "Дополнительное сечение рода 1" M=m треугольников=polygons пикселей=count(mask) #???
        # !!! keep_masks && push!(masks, mask)
        push!(masks, mask) #???
        combined .|= mask #???
    end #???
    return masks, combined, τ #???
end #???

# !!! Для преобразованного графика каждая рыжая часть — один треугольник.
# !!! Для преобразованного графика каждая зелёная часть — один треугольник.
#??? Для преобразованного графика каждая средне-серая часть — один треугольник
#??? с вершинами на диагонали (a,a), (b,b) и в углу (a,b).
function _slice_gen1_triangles_transformed(z0; N, mobius, coordinates) #???
    τ = _slice_gen1_tau(z0, N) #???
    right = Float64[x_gen_1(τ, N; right_point=k) for k in 0:N] #???
    transformed = [mobius.transform(value) for value in right] #???
    pixels = length(coordinates) #???
    triangles = Vector{Vector{_SlicePoint}}() #???
    masks = BitMatrix[] #???
    combined = falses(pixels, pixels) #???
    for k in 1:N #???
        a, b = transformed[k], transformed[k + 1] #???
        if !isfinite(a) || !isfinite(b) || b <= a #???
            # !!! @warn "Рыжий треугольник пропущен: неверный порядок граничных точек" k=k a=a b=b
            # !!! @warn "Зелёный треугольник пропущен: неверный порядок граничных точек" k=k a=a b=b
            @warn "Треугольник рода 1 пропущен: неверный порядок граничных точек" k=k a=a b=b #???
            continue #???
        end #???
        a, b = clamp(a, -1.0, 1.0), clamp(b, -1.0, 1.0) #???
        b > a || continue #???
        triangle = _SlicePoint[(a, a), (a, b), (b, b)] #???
        mask = falses(pixels, pixels) #???
        _slice_paint!(mask, coordinates, triangle) #???
        push!(triangles, triangle) #???
        push!(masks, mask) #???
        combined .|= mask #???
    end #???
    # !!! @info "Рыжие треугольники по трём вершинам" всего=length(triangles) ожидается=N пикселей=count(combined)
    # !!! @info "Зелёные треугольники по трём вершинам" всего=length(triangles) ожидается=N пикселей=count(combined)
    @info "Треугольники рода 1 по трём вершинам" всего=length(triangles) ожидается=N пикселей=count(combined) #???
    return triangles, masks, combined, τ #???
end #???

# !!! Рисуем рыжие области непосредственно по трём вершинам, без пиксельной заливки.
# !!! Рисуем зелёные области непосредственно по трём вершинам, без пиксельной заливки.
#??? Рисуем средне-серые области непосредственно по трём вершинам.
function _slice_draw_gen1_triangles!(axis, triangles; strokewidth=2.0) #???
    # !!! fill_color = Makie.RGBAf(1.0, 0.50, 0.08, 1.0)
    # !!! edge_color = Makie.RGBAf(0.7, 0.22, 0.0, 1.0)
    # !!! fill_color = Makie.RGBAf(0.20, 0.72, 0.32, 1.0)
    # !!! edge_color = Makie.RGBAf(0.08, 0.43, 0.16, 1.0)
    fill_color = Makie.RGBAf(0.62, 0.62, 0.62, 1.0) #???
    edge_color = Makie.RGBAf(0.27, 0.27, 0.27, 1.0) #???
    for triangle in triangles #???
        points = [Makie.Point2f(x, y) for (x, y) in triangle] #???
        Makie.poly!(axis, points; color=fill_color, #???
                    strokecolor=edge_color, strokewidth=strokewidth) #???
    end #???
    return nothing #???
end #???

#= !!! Зелёная основа для незакрашенного остатка области X<Y после преобразования.
function _slice_draw_green_background!(axis) #???
    triangle = Makie.Point2f[Makie.Point2f(-1.0, -1.0), #???
                             Makie.Point2f(-1.0, 1.0), #???
                             Makie.Point2f(1.0, 1.0)] #???
    Makie.poly!(axis, triangle; color=Makie.RGBAf(0.20, 0.72, 0.32, 1.0)) #???
    return nothing #???
end #???
=#
#= !!! Рыжая основа для остатка области X<Y; голубая и зелёная заливки поверх неё.
function _slice_draw_orange_background!(axis) #???
    triangle = Makie.Point2f[Makie.Point2f(-1.0, -1.0), #???
                             Makie.Point2f(-1.0, 1.0), #???
                             Makie.Point2f(1.0, 1.0)] #???
    Makie.poly!(axis, triangle; color=Makie.RGBAf(1.0, 0.50, 0.08, 1.0)) #???
    return nothing #???
end #???
=#
#??? Светло-серая основа для остатка области X<Y; другие области поверх неё.
function _slice_draw_gray_background!(axis) #???
    triangle = Makie.Point2f[Makie.Point2f(-1.0, -1.0), #???
                             Makie.Point2f(-1.0, 1.0), #???
                             Makie.Point2f(1.0, 1.0)] #???
    Makie.poly!(axis, triangle; color=Makie.RGBAf(0.86, 0.86, 0.86, 1.0)) #???
    return nothing #???
end #???

# !!! Рисует сечение двух областей плоскостью `Z=z0`. Все области закрашены одним цветом.
# !!! Тонкие линии показывают границы отдельных `M`, толстая линия — общую границу.
#??? Описание дополнено дополнительной областью рода 1.
"""
    plot_volume_slice(z0; N=8, M=1:(N-1), x_bounds=(1e-3, 10.0),
                        y_bounds=(1e-3, 10.0), nx=31, ny=31, pixels=500,
                        refine_intersections=true)

Рисует сечение прежних областей и двух треугольников рода 1 для каждого `M`.
Тонкие линии показывают границы отдельных `M`, толстые — общие границы.
# !!! `refine_intersections=true` уточняло пересечения x₆=z0 методом Ньютона с ограничением внутри ребра.
# !!! `refine_intersections=true` уточняло пересечения методом секущих внутри ребра.
`refine_intersections=true` уточняет пересечения x₆=z0 методом Ньютона внутри ребра; при неудаче используется линейная оценка. #???
"""
# !!! function plot_volume_slice(z0; N=8, M=1:(N-1), x_bounds=(1e-3, 10.0), y_bounds=(1e-3, 10.0), nx=31, ny=31, pixels=500)
function plot_volume_slice(z0; N=8, M=1:(N-1), x_bounds=(1e-3, 10.0),
                           y_bounds=(1e-3, 10.0), nx=31, ny=31, pixels=500,
                           refine_intersections=true) #???
    # !!! 0 < z0 < 100 || throw(ArgumentError("нужно 0 < z0 < 100"))
    # !!! 0 < z0 < 100 || throw(ArgumentError("нужно 0 < z0 < 100; получено z0=$(repr(z0))"))
    isfinite(z0) && z0 > 0 || throw(ArgumentError("нужно конечное z0 > 0; получено z0=$(repr(z0))")) #???
    N isa Integer && N >= 2 || throw(ArgumentError("N должно быть целым числом не меньше 2"))
    0 < x_bounds[1] < x_bounds[2] || throw(ArgumentError("нужно 0 < x_bounds[1] < x_bounds[2]"))
    0 < y_bounds[1] < y_bounds[2] || throw(ArgumentError("нужно 0 < y_bounds[1] < y_bounds[2]"))
    nx >= 2 && ny >= 2 || throw(ArgumentError("nx и ny должны быть не меньше 2"))
    pixels >= 2 || throw(ArgumentError("pixels должно быть не меньше 2"))
    M_values = M isa Integer ? [M] : collect(M)
    isempty(M_values) && throw(ArgumentError("список M не может быть пустым"))
    all(m -> m isa Integer && 1 <= m < N, M_values) ||
        throw(ArgumentError("каждое M должно быть целым числом из 1:N-1"))
    M_values = unique(Int.(M_values))

    masks = BitMatrix[]
    combined = falses(pixels, pixels)
    for m in M_values
        # !!! mask, polygons = _slice_mask_for_M(z0; N=N, M=m, x_bounds=x_bounds, y_bounds=y_bounds, nx=nx, ny=ny, pixels=pixels)
        mask, polygons = _slice_mask_for_M(z0; N=N, M=m, x_bounds=x_bounds, #???
                                          y_bounds=y_bounds, nx=nx, ny=ny, pixels=pixels, #???
                                          refine_intersections=refine_intersections) #???
        @info "Сечение" M=m многоугольников=polygons пикселей=count(mask)
        push!(masks, mask)
        combined .|= mask
    end
    gen1_masks, gen1_combined, τ = _slice_gen1_masks(z0; N=N, M_values=M_values, pixels=pixels) #???
    @info "Высота прямоугольника рода 1 для сечения" Z=z0 τ=τ #???
    # !!! count(combined) == 0 && @warn "Сечение пусто; измените z0 или диапазоны s,t."
    if count(combined) == 0 && count(gen1_combined) == 0 #???
        @warn "Сечение пусто; измените z0 или диапазоны s,t." #???
    end #???

    coordinates = range(0.0, Float64(z0); length=pixels)
    fig = Makie.Figure(size=(850, 800))
    # !!! ax = Makie.Axis(fig[1, 1]; xlabel="X", ylabel="Y", aspect=Makie.DataAspect(),
    # !!!                 title="Сечение Z = $(z0)")
    ax = Makie.Axis(fig[1, 1]; xlabel="X", ylabel="Y", aspect=Makie.DataAspect(), #???
                    title="Z = $(z0): синий — прежние области, оранжевый — род 1") #???
    fill_color = Makie.RGBAf(0.12, 0.56, 1.0, 1.0)
    Makie.heatmap!(ax, coordinates, coordinates, Float32.(combined);
                   colormap=[Makie.RGBAf(1, 1, 1, 0), fill_color], colorrange=(0, 1))
    if any(gen1_combined) #???
        gen1_color = Makie.RGBAf(1.0, 0.50, 0.08, 0.78) #???
        Makie.heatmap!(ax, coordinates, coordinates, Float32.(gen1_combined); #???
                       colormap=[Makie.RGBAf(1, 1, 1, 0), gen1_color], colorrange=(0, 1)) #???
    end #???
    for mask in masks
        any(mask) || continue
        Makie.contour!(ax, coordinates, coordinates, Float32.(mask);
                       levels=[0.5], color=(:black, 0.35), linewidth=1)
    end
    if any(combined)
        Makie.contour!(ax, coordinates, coordinates, Float32.(combined);
                       levels=[0.5], color=:black, linewidth=2.5)
    end
    for mask in gen1_masks #???
        any(mask) || continue #???
        Makie.contour!(ax, coordinates, coordinates, Float32.(mask); #???
                       levels=[0.5], color=Makie.RGBAf(0.7, 0.22, 0.0, 0.5), linewidth=1) #???
    end #???
    if any(gen1_combined) #???
        Makie.contour!(ax, coordinates, coordinates, Float32.(gen1_combined); #???
                       levels=[0.5], color=Makie.RGBAf(0.7, 0.22, 0.0, 1.0), linewidth=2.5) #???
    end #???
    Makie.lines!(ax, [0.0, 0.0, z0, 0.0], [0.0, z0, z0, 0.0];
                 color=(:gray, 0.6), linestyle=:dash)
    Makie.xlims!(ax, 0, z0)
    Makie.ylims!(ax, 0, z0)
    display(fig)
    return fig
end

#??? Устойчивая запись преобразования t(s)=(s+a)/(Ls-a).
#??? Здесь L=z/(sqrt(1+z)+1)^2 и a=-z/(sqrt(1+z)+1) равны формулам через 1+2/z.
function slice_mobius_transform(z) #???
    isfinite(z) && z > 0 || throw(ArgumentError("нужно конечное z > 0")) #???
    root = sqrt(1 + Float64(z)) #???
    L = Float64(z) / (root + 1)^2 #???
    a = -Float64(z) / (root + 1) #???
    transform(s) = isinf(s) ? inv(L) : (s + a) / (L * s - a) #???
    inverse(u) = a * (u + 1) / (L * u - 1) #???
    return (; transform, inverse, L, a) #???
end #???

#= !!! Отражение (u, v) -> (-v, -u) на квадратной сетке: (i, j) -> (n+1-j, n+1-i).
function _slice_symmetrize_upper_right(mask::AbstractMatrix{Bool}) #???
    n, m = size(mask) #???
    n == m || throw(ArgumentError("маска должна быть квадратной")) #???
    symmetric_mask = copy(mask) #???
    mismatched_pairs = 0 #???
    for j in 1:n, i in 1:n #???
        if i + j < n + 1 #???
            reflected = mask[n + 1 - j, n + 1 - i] #???
            mismatched_pairs += mask[i, j] != reflected #???
            symmetric_mask[i, j] = reflected #???
        end #???
    end #???
    return symmetric_mask, mismatched_pairs #???
end #???
=#

#??? Отражаем верхнюю правую половину относительно t(Y)=-t(X): (u,v) ↦ (-v,-u).
function _slice_reflect_upper_right(mask::AbstractMatrix{Bool}) #???
    n, m = size(mask) #???
    n == m || throw(ArgumentError("маска должна быть квадратной")) #???
    result = copy(mask) #???
    for j in 1:n, i in 1:n #???
        if i + j < n + 1 #???
            result[i, j] = mask[n + 1 - j, n + 1 - i] #???
        end #???
    end #???
    return result #???
end #???

# !!! Нижняя левая половина строится отражением верхней правой относительно
# !!! `t(Y)=-t(X)`; число несовпадающих исходных пикселей выводится в журнал.
# !!! Прежняя документация: M=1:(N-1), обе половины рассчитывались непосредственно.
# !!! Прежний отражённый график рассчитывал M=N-1, …, ceil(N/2).
# !!! Голубая нижняя часть получалась отражением верхней относительно t(Y)=-t(X).
""" #???
    plot_volume_slice_transformed(z0; N=8, M=nothing,
                                  x_bounds=(1e-3, 10.0), y_bounds=(1e-3, 10.0),
                                  nx=31, ny=31, pixels=500,
                                  refine_intersections=true,
                                  individual_boundaries=false,
                                  lower_zoom=true, zoom_limit=-0.5,
                                  graph_name="volume_slice_transformed") #???

Рисует то же сечение Z=z0 после преобразования каждой координаты
`t(s)=(s+a)/(L*s-a)`. Новые координаты удовлетворяют `-1 < t(X) < t(Y) < 1`.
Строит маски на равномерной сетке преобразованных координат,
чтобы мелкие области у края t=-1 оставались различимыми.
`refine_intersections=true` уточняет пересечения x₆=z0 на рёбрах параметрической сетки с линейной оценкой при неудаче. #???
Границы отдельных M можно включить через `individual_boundaries=true`.
# !!! `lower_zoom=true` добавляет увеличенный вид области около (-1,-1).
# !!! Рисуется только увеличенный вид области около (-1,-1); `zoom_limit` задаёт его верхнюю границу.
Рисуется только общий вид прежней левой панели на квадрате `[-1,1]²`. #???
`lower_zoom` и `zoom_limit` сохранены для совместимости со старыми вызовами. #???
`graph_name` задаёт имя PNG-файла рядом с `volume_slice.jl`; расширение `.png` можно опустить. #???
# !!! Голубые области рассчитываются при `M=1, …, floor(N/2)` без отражения.
# !!! Голубые области рассчитываются при `M=1, …, floor(N/2)`; верхняя правая часть отражается относительно `t(Y)=-t(X)`.
Тёмно-серые области рассчитываются при `M=1, …, floor(N/2)`; верхняя правая часть отражается относительно `t(Y)=-t(X)`. #???
# !!! Рыжая область рода 1 рассчитывается для всех `M=1, …, N-1` без отражения.
# !!! Зелёная область рода 1 рассчитывается для всех `M=1, …, N-1` без отражения.
Средне-серые треугольники рода 1 рассчитываются для всех `M=1, …, N-1` без отражения. #???
# !!! Остальная часть треугольника `t(X)<t(Y)` закрашивается зелёным.
# !!! Остальная часть треугольника `t(X)<t(Y)` закрашивается рыжим.
#??? Остальная часть треугольника `t(X)<t(Y)` закрашивается светло-серым.
#??? Каждый её участок рисуется одним треугольником по трём вершинам в координатах t.
Точки сетки лежат строго внутри `-1 < t(X), t(Y) < 1`. #???
""" #???
# !!! function plot_volume_slice_transformed(z0; N=8, M=nothing, x_bounds=(1e-3, 10.0), y_bounds=(1e-3, 10.0), nx=31, ny=31, pixels=500, individual_boundaries=false, lower_zoom=true, zoom_limit=-0.5)
function plot_volume_slice_transformed(z0; N=8, M=nothing, #???
                                       x_bounds=(1e-3, 10.0), y_bounds=(1e-3, 10.0), #???
                                       nx=31, ny=31, pixels=500,
                                       refine_intersections=true, #???
                                       individual_boundaries=false,
                                       # !!! lower_zoom=true, zoom_limit=-0.5)
                                       lower_zoom=true, zoom_limit=-0.5,
                                       graph_name::AbstractString="volume_slice_transformed") #???
    mobius = slice_mobius_transform(z0) #???
    N isa Integer && N >= 2 || throw(ArgumentError("N должно быть целым числом не меньше 2")) #???
    0 < x_bounds[1] < x_bounds[2] || throw(ArgumentError("нужно 0 < x_bounds[1] < x_bounds[2]")) #???
    0 < y_bounds[1] < y_bounds[2] || throw(ArgumentError("нужно 0 < y_bounds[1] < y_bounds[2]")) #???
    nx >= 2 && ny >= 2 || throw(ArgumentError("nx и ny должны быть не меньше 2")) #???
    pixels >= 2 || throw(ArgumentError("pixels должно быть не меньше 2")) #???
    -1 < zoom_limit < 1 || throw(ArgumentError("нужно -1 < zoom_limit < 1")) #???
    isempty(strip(graph_name)) && throw(ArgumentError("graph_name не может быть пустым")) #???
    # !!! M_values = M isa Integer ? [M] : collect(M) #???
    # !!! M_values = isnothing(M) ? collect((N - 1):-1:cld(N, 2)) : ...
    M_values = isnothing(M) ? collect(1:fld(N, 2)) : #???
               (M isa Integer ? [M] : collect(M)) #???
    isempty(M_values) && throw(ArgumentError("список M не может быть пустым")) #???
    all(m -> m isa Integer && 1 <= m < N, M_values) || #???
        throw(ArgumentError("каждое M должно быть целым числом из 1:N-1")) #???
    # !!! M_values = unique(Int.(M_values)) #???
    # !!! M_values = sort!(unique(Int.(M_values)); rev=true) #???
    M_values = sort!(unique(Int.(M_values))) #???
    # !!! all(m -> m >= cld(N, 2), M_values) || throw(ArgumentError("для отражённого графика нужно M ≥ ceil(N/2)")) #???
    all(m -> m <= fld(N, 2), M_values) || #???
        throw(ArgumentError("для этого графика нужно M ≤ floor(N/2)")) #???
    # !!! @info "Расчёт верхней половины" M=M_values #???
    # !!! @info "Расчёт голубых областей без отражения" M=M_values
    # !!! @info "Расчёт голубых областей с отражением" M=M_values
    @info "Расчёт основных областей с отражением" M=M_values #???

    # !!! transformed_coordinates = range(-1.0, 1.0; length=pixels) #???
    #??? Центры ячеек лежат строго внутри (-1,1); сами t=-1 и t=1 не вычисляются.
    transformed_coordinates = range(-1.0 + 1.0 / pixels, 1.0 - 1.0 / pixels; length=pixels) #???
    #= !!! Прежний способ переносил уже готовую маску с равномерной сетки 0...z0.
    # !!! Из-за повторения исходных пикселей мелкие области около t=-1 становились ступенчатыми.
    # !!! original_coordinates = [mobius.inverse(u) for u in transformed_coordinates] #???
    # !!! original_coordinates[1] = 0.0 #???
    # !!! original_coordinates[end] = Float64(z0) #???
    # !!! masks = BitMatrix[] #???
    # !!! combined = falses(pixels, pixels) #???
    original_coordinates = range(0.0, Float64(z0); length=pixels) #???
    source_masks = BitMatrix[] #???
    source_combined = falses(pixels, pixels) #???
    for m in M_values #???
        # !!! mask, polygons = _slice_mask_for_M(z0; N=N, M=m, x_bounds=x_bounds, y_bounds=y_bounds, nx=nx, ny=ny, pixels=pixels, coordinates=original_coordinates)
        mask, polygons = _slice_mask_for_M(z0; N=N, M=m, x_bounds=x_bounds, #???
                                          y_bounds=y_bounds, nx=nx, ny=ny, #???
                                          pixels=pixels, coordinates=original_coordinates, #???
                                          refine_intersections=refine_intersections) #???
        @info "Исходное сечение для преобразования" M=m многоугольников=polygons пикселей=count(mask) #???
        # !!! push!(masks, mask) #???
        # !!! combined .|= mask #???
        push!(source_masks, mask) #???
        source_combined .|= mask #???
    end #???
    # !!! gen1_masks, gen1_combined, τ = _slice_gen1_masks(z0; N=N, M_values=M_values, #???
    # !!!                                              pixels=pixels, coordinates=original_coordinates) #???
    source_gen1_masks, source_gen1_combined, τ = _slice_gen1_masks( #???
        z0; N=N, M_values=M_values, pixels=pixels, coordinates=original_coordinates) #???
    #??? Берём те же пиксели, что и на исходном графике, затем меняем только их координаты.
    preimage_coordinates = [mobius.inverse(u) for u in transformed_coordinates] #???
    source_indices = [clamp(round(Int, s / z0 * (pixels - 1) + 1), 1, pixels) #???
                      for s in preimage_coordinates] #???
    # !!! masks = [mask[source_indices, source_indices] for mask in source_masks] #???
    masks = individual_boundaries ? [mask[source_indices, source_indices] for mask in source_masks] : BitMatrix[] #???
    combined = source_combined[source_indices, source_indices] #???
    gen1_masks = [mask[source_indices, source_indices] for mask in source_gen1_masks] #???
    gen1_combined = source_gen1_combined[source_indices, source_indices] #???
    =#
    #??? Каждый пиксель нового графика проверяем в его точном прообразе, без округления к старому пикселю.
    original_coordinates = [mobius.inverse(u) for u in transformed_coordinates] #???
    # !!! original_coordinates[1] = 0.0 #???
    # !!! original_coordinates[end] = Float64(z0) #???
    masks = BitMatrix[] #???
    combined = falses(pixels, pixels) #???
    for m in M_values #???
        # !!! mask, polygons = _slice_mask_for_M(z0; N=N, M=m, x_bounds=x_bounds, y_bounds=y_bounds, nx=nx, ny=ny, pixels=pixels, coordinates=original_coordinates)
        mask, polygons = _slice_mask_for_M(z0; N=N, M=m, x_bounds=x_bounds, #???
                                          y_bounds=y_bounds, nx=nx, ny=ny, #???
                                          pixels=pixels, coordinates=original_coordinates, #???
                                          refine_intersections=refine_intersections) #???
        # !!! @info "Преобразованное сечение" M=m многоугольников=polygons пикселей=count(mask) #???
        upper_pixels = count(mask[i, j] && i + j >= pixels + 1 #???
                             for j in 1:pixels, i in 1:pixels) #???
        @info "Преобразованное сечение" M=m многоугольников=polygons пикселей=count(mask) над_диагональю=upper_pixels #???
        push!(masks, mask) #???
        combined .|= mask #???
    end #???
    #= !!! Рыжая область ошибочно ограничивалась голубым диапазоном M.
    gen1_masks, gen1_combined, τ = _slice_gen1_masks(
        z0; N=N, M_values=M_values, pixels=pixels, coordinates=original_coordinates)
    =#
    #= !!! Рыжие области заполнялись по маске из двух треугольников для каждого M.
    gen1_M_values = collect(1:(N - 1))
    gen1_masks, gen1_combined, τ = _slice_gen1_masks(
        z0; N=N, M_values=gen1_M_values, pixels=pixels, coordinates=original_coordinates)
    =#
    gen1_triangles, gen1_masks, gen1_combined, τ = _slice_gen1_triangles_transformed( #???
        z0; N=N, mobius=mobius, coordinates=transformed_coordinates) #???
    # !!! Отражаем только голубые маски; рыжие треугольники уже построены для всех интервалов.
    # !!! Отражаем только голубые маски; зелёные треугольники уже построены для всех интервалов.
    #??? Отражаем только маски основных областей; треугольники рода 1 не отражаем.
    masks = [_slice_reflect_upper_right(mask) for mask in masks] #???
    combined = _slice_reflect_upper_right(combined) #???
    # !!! Голубые маски остаются в исходном положении; рыжие треугольники тоже не отражаются.
    # !!! gen1_masks = [_slice_reflect_upper_right(mask) for mask in gen1_masks] #???
    # !!! gen1_combined = _slice_reflect_upper_right(gen1_combined) #???
    # !!! @info "Совпадение синих и оранжевых масок" пикселей=count(combined .& gen1_combined)
    # !!! @info "Совпадение синих и зелёных масок" пикселей=count(combined .& gen1_combined)
    @info "Совпадение основных масок и масок рода 1" пикселей=count(combined .& gen1_combined) #???
    #??? Диагностика совпадений относится к текущей, равномерной по t сетке.
    occupied = falses(pixels, pixels) #???
    source_overlap = falses(pixels, pixels) #???
    # !!! for mask in source_masks #???
    for mask in masks #???
        source_overlap .|= occupied .& mask #???
        occupied .|= mask #???
    end #???
    # !!! @info "Пересечение исходных синих масок" пикселей=count(source_overlap) #???
    # !!! if any(source_overlap) #???
    # !!!     @warn "Синие маски для разных M пересекаются уже на исходном сечении; преобразование сохраняет это пересечение." #???
    # !!! end #???
    #??? Одинокие общие пиксели на границе отличаются от перекрытия с ненулевой толщиной.
    overlap_interior = falses(pixels, pixels) #???
    if pixels >= 3 #???
        overlap_interior[2:end-1, 2:end-1] .= #???
            source_overlap[2:end-1, 2:end-1] .& #???
            source_overlap[1:end-2, 2:end-1] .& #???
            source_overlap[3:end, 2:end-1] .& #???
            source_overlap[2:end-1, 1:end-2] .& #???
            source_overlap[2:end-1, 3:end] #???
    end #???
    # !!! @info "Совпадение синих масок" общих_пикселей=count(source_overlap)
    @info "Совпадение основных масок" общих_пикселей=count(source_overlap) #???
    @info "Совпадение вне однопиксельной границы" внутренних_пикселей=count(overlap_interior) #???
    if any(overlap_interior) #???
        # !!! @warn "Области для разных M перекрываются на исходном сечении вне однопиксельной границы." #???
        @warn "Области для разных M перекрываются на преобразованной сетке вне однопиксельной границы." #???
    end #???
    #= !!! Прежнее отражение всех рассчитанных M, без ограничения верхней половины.
    combined, blue_mismatch = _slice_symmetrize_upper_right(combined) #???
    gen1_combined, orange_mismatch = _slice_symmetrize_upper_right(gen1_combined) #???
    @info "Асимметрия исходного сечения относительно t(Y)=-t(X)" синих_пар=blue_mismatch оранжевых_пар=orange_mismatch #???
    masks = [_slice_symmetrize_upper_right(mask)[1] for mask in masks] #???
    =#
    # !!! Ранее отражение выполнялось после диагностики пересечений.
    # !!! Оранжевый имеет приоритет там, где исходные маски совпадают.
    # !!! Зелёные треугольники имеют приоритет там, где маски совпадают.
    #??? Треугольники рода 1 имеют приоритет там, где маски совпадают.
    # !!! blue_before_orange = count(combined)
    # !!! blue_before_green = count(combined)
    main_before_gen1 = count(combined) #???
    combined .&= .!gen1_combined #???
    masks = [mask .& .!gen1_combined for mask in masks] #???
    # !!! @info "Голубая область после наложения оранжевой" до=blue_before_orange после=count(combined)
    # !!! @info "Голубая область после наложения зелёной" до=blue_before_green после=count(combined)
    @info "Основная область после наложения рода 1" до=main_before_gen1 после=count(combined) #???
    @info "Высота прямоугольника рода 1 для преобразованного сечения" Z=z0 τ=τ #???
    if !any(combined) && !any(gen1_combined) #???
        @warn "Сечение пусто; измените z0 или диапазоны s,t." #???
    end #???

    #= !!! Прежний рисунок с общим видом слева и увеличением справа.
    # !!! fig = Makie.Figure(size=(850, 800)) #???
    fig = Makie.Figure(size=lower_zoom ? (1650, 850) : (850, 800)) #???
    # !!! ax = Makie.Axis(fig[1, 1]; xlabel="t(X)", ylabel="t(Y)", aspect=Makie.DataAspect(), #???
    # !!!                 title="t(Z) = 1 (Z = $(z0)): синий — прежние области, оранжевый — род 1") #???
    ax = Makie.Axis(fig[1, 1]; xlabel="t(X)", ylabel="t(Y)", aspect=Makie.DataAspect(), #???
                    title="Сечение t(Z) = 1 (Z = $(z0))") #???
    # !!! _slice_draw_green_background!(ax)
    # !!! _slice_draw_orange_background!(ax)
    _slice_draw_gray_background!(ax) #???
    # !!! fill_color = Makie.RGBAf(0.12, 0.56, 1.0, 1.0)
    fill_color = Makie.RGBAf(0.30, 0.30, 0.30, 1.0) #???
    Makie.heatmap!(ax, transformed_coordinates, transformed_coordinates, Float32.(combined); #???
                   colormap=[Makie.RGBAf(1, 1, 1, 0), fill_color], colorrange=(0, 1)) #???
    #= !!! Рыжая заливка строилась через пиксельную маску.
    if any(gen1_combined) #???
        # !!! gen1_color = Makie.RGBAf(1.0, 0.50, 0.08, 0.78) #???
        gen1_color = Makie.RGBAf(1.0, 0.50, 0.08, 1.0) #???
        Makie.heatmap!(ax, transformed_coordinates, transformed_coordinates, Float32.(gen1_combined); #???
                       colormap=[Makie.RGBAf(1, 1, 1, 0), gen1_color], colorrange=(0, 1)) #???
    end #???
    =#
    # !!! for mask in masks #???
    # !!!     any(mask) || continue #???
    # !!!     Makie.contour!(ax, transformed_coordinates, transformed_coordinates, Float32.(mask); #???
    # !!!                    levels=[0.5], color=(:black, 0.35), linewidth=1) #???
    # !!! end #???
    if individual_boundaries #???
        for mask in masks #???
            any(mask) || continue #???
            Makie.contour!(ax, transformed_coordinates, transformed_coordinates, Float32.(mask); #???
                           levels=[0.5], color=(:black, 0.35), linewidth=1) #???
        end #???
    end #???
    if any(combined) #???
        Makie.contour!(ax, transformed_coordinates, transformed_coordinates, Float32.(combined); #???
                       levels=[0.5], color=:black, linewidth=2.5) #???
    end #???
    _slice_draw_gen1_triangles!(ax, gen1_triangles; strokewidth=2.5) #???
    #= !!! Внутренние границы оранжевых масок разных M затемняли заливку.
    for mask in gen1_masks #???
        any(mask) || continue #???
        Makie.contour!(ax, transformed_coordinates, transformed_coordinates, Float32.(mask); #???
                       levels=[0.5], color=Makie.RGBAf(0.7, 0.22, 0.0, 0.5), linewidth=1) #???
    end #???
    =#
    #= !!! Рыжая граница повторяла ступени пиксельной маски.
    if any(gen1_combined) #???
        Makie.contour!(ax, transformed_coordinates, transformed_coordinates, Float32.(gen1_combined); #???
                       levels=[0.5], color=Makie.RGBAf(0.7, 0.22, 0.0, 1.0), linewidth=2.5) #???
    end #???
    =#
    Makie.lines!(ax, [-1.0, -1.0, 1.0, -1.0], [-1.0, 1.0, 1.0, -1.0]; #???
                 color=(:gray, 0.6), linestyle=:dash) #???
    Makie.xlims!(ax, -1, 1) #???
    Makie.ylims!(ax, -1, 1) #???
    if lower_zoom #???
        detail = Makie.Axis(fig[1, 2]; xlabel="t(X)", ylabel="t(Y)", #???
                            aspect=Makie.DataAspect(), title="Вблизи (-1, -1)") #???
        # !!! _slice_draw_green_background!(detail)
        # !!! _slice_draw_orange_background!(detail)
        _slice_draw_gray_background!(detail) #???
        Makie.heatmap!(detail, transformed_coordinates, transformed_coordinates, Float32.(combined); #???
                       colormap=[Makie.RGBAf(1, 1, 1, 0), fill_color], colorrange=(0, 1)) #???
        #= !!! Рыжая заливка на увеличении строилась через пиксельную маску.
        if any(gen1_combined) #???
            # !!! gen1_color = Makie.RGBAf(1.0, 0.50, 0.08, 0.78) #???
            gen1_color = Makie.RGBAf(1.0, 0.50, 0.08, 1.0) #???
            Makie.heatmap!(detail, transformed_coordinates, transformed_coordinates, Float32.(gen1_combined); #???
                           colormap=[Makie.RGBAf(1, 1, 1, 0), gen1_color], colorrange=(0, 1)) #???
        end #???
        =#
        if any(combined) #???
            Makie.contour!(detail, transformed_coordinates, transformed_coordinates, Float32.(combined); #???
                           levels=[0.5], color=:black, linewidth=2) #???
        end #???
        _slice_draw_gen1_triangles!(detail, gen1_triangles; strokewidth=2.0) #???
        #= !!! Рыжая граница на увеличении повторяла ступени пиксельной маски.
        if any(gen1_combined) #???
            Makie.contour!(detail, transformed_coordinates, transformed_coordinates, Float32.(gen1_combined); #???
                           levels=[0.5], color=Makie.RGBAf(0.7, 0.22, 0.0, 1.0), linewidth=2) #???
        end #???
        =#
        Makie.lines!(detail, [-1.0, -1.0, zoom_limit, -1.0], #???
                     [-1.0, zoom_limit, zoom_limit, -1.0]; #???
                     color=(:gray, 0.6), linestyle=:dash) #???
        Makie.xlims!(detail, -1, zoom_limit) #???
        Makie.ylims!(detail, -1, zoom_limit) #???
    end #???
    display(fig) #???
    return fig #???
    =#

    # !!! fig = Makie.Figure(size=(850, 800)) #??? Оставляем только прежнюю правую панель.
    # !!! fig = Makie.Figure(size=(850, 800)) #??? Оставляем только прежнюю левую панель.
    fig = Makie.Figure(size=(1400, 1300)) #??? Увеличенный рисунок прежней левой панели.
    # !!! Makie.Label(fig[0, 1], "x(p₆) = $(z0)"; fontsize=24)
    Makie.Label(fig[0, 1], "λ = $(z0)"; fontsize=32) #??? Греческая λ и уровень Z=z0.
    ax = Makie.Axis(fig[1, 1]; aspect=Makie.DataAspect()) #??? Без подписей осей и старых заголовков.
    Makie.colsize!(fig.layout, 1, Makie.Relative(1)) #??? Ось занимает ширину рисунка, а не ширину подписи.
    Makie.rowsize!(fig.layout, 1, Makie.Relative(1)) #??? Ось занимает высоту под верхней подписью.
    _slice_draw_gray_background!(ax) #???
    fill_color = Makie.RGBAf(0.30, 0.30, 0.30, 1.0) #???
    Makie.heatmap!(ax, transformed_coordinates, transformed_coordinates, Float32.(combined); #???
                   colormap=[Makie.RGBAf(1, 1, 1, 0), fill_color], colorrange=(0, 1)) #???
    if individual_boundaries #???
        for mask in masks #???
            any(mask) || continue #???
            Makie.contour!(ax, transformed_coordinates, transformed_coordinates, Float32.(mask); #???
                           levels=[0.5], color=(:black, 0.35), linewidth=1) #???
        end #???
    end #???
    if any(combined) #???
        # !!! Makie.contour!(ax, transformed_coordinates, transformed_coordinates, Float32.(combined); levels=[0.5], color=:black, linewidth=2)
        Makie.contour!(ax, transformed_coordinates, transformed_coordinates, Float32.(combined); #???
                       levels=[0.5], color=:black, linewidth=2.5) #???
    end #???
    # !!! _slice_draw_gen1_triangles!(ax, gen1_triangles; strokewidth=2.0)
    _slice_draw_gen1_triangles!(ax, gen1_triangles; strokewidth=2.5) #???
    # !!! Makie.lines!(ax, [-1.0, -1.0, zoom_limit, -1.0], [-1.0, zoom_limit, zoom_limit, -1.0]; color=(:gray, 0.6), linestyle=:dash)
    Makie.lines!(ax, [-1.0, -1.0, 1.0, -1.0], #???
                 [-1.0, 1.0, 1.0, -1.0]; #???
                 color=(:gray, 0.6), linestyle=:dash) #???
    # !!! Makie.xlims!(ax, -1, zoom_limit)
    # !!! Makie.ylims!(ax, -1, zoom_limit)
    Makie.xlims!(ax, -1, 1) #???
    Makie.ylims!(ax, -1, 1) #???
    Makie.hidedecorations!(ax) #??? Убираем подписи осей, делений и сетку.
    graph_file = endswith(lowercase(graph_name), ".png") ? graph_name : graph_name * ".png" #???
    graph_path = joinpath(@__DIR__, graph_file) #??? Сохраняем рядом с volume_slice.jl.
    Makie.save(graph_path, fig) #???
    println("График сохранён: ", graph_path) #???
    display(fig) #???
    return fig #???
end #???

