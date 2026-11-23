#??? Новый независимый график: сечение при фиксированном γ = T(x₁).
#??? При прямом include подключаем определения xᵢ из points.jl.
if !isdefined(@__MODULE__, :x_1) #???
    include(joinpath(@__DIR__, "points.jl")) #???
end #???

const _GammaPoint = NTuple{2, Float64} #???

#= !!! Нули многочлена Чебышёва в логарифмической координате вместе с
#!!! концами интервала. Прежняя сетка с концентрацией у краёв.
function _gamma_chebyshev_log_nodes(bounds, n::Int) #???
    n >= 2 || throw(ArgumentError("нужно хотя бы два узла")) #???
    low, high = log(Float64(bounds[1])), log(Float64(bounds[2])) #???
    nodes = Vector{Float64}(undef, n) #???
    nodes[1], nodes[end] = Float64(bounds[1]), Float64(bounds[2]) #???
    roots = n - 2 #???
    for k in 1:roots #???
        θ = (2k - 1) * π / (2roots) #???
        nodes[k + 1] = exp((low + high) / 2 - #???
                           (high - low) * cos(θ) / 2) #???
    end #???
    return nodes #???
end #???
=#

function _gamma_finite_value(f) #???
    try #???
        value = f() #???
        return isfinite(value) ? Float64(value) : NaN #???
    catch err #???
        if err isa DomainError || err isa DivideError || err isa OverflowError || #???
           err isa ArgumentError || err isa LinearAlgebra.SingularException || #???
           err isa LinearAlgebra.PosDefException #???
            return NaN #???
        end #???
        rethrow() #???
    end #???
end #???

#??? T(x₃)=-1, T(x₆)=1, T(x₁)=-T(x₂).
# !!! Раньше запрещали полюс между x₃ и x₆ и выбирали неверную ветвь.
#??? Из двух вещественных ветвей выбираем ту, где -1<T(x₄)<T(x₅)<1.
function _gamma_normalization(points::NTuple{6, Float64}) #???
    all(isfinite, points) && length(unique(points)) == 6 || return nothing #???
    x1, x2, x3, x4, x5, x6 = points #???
    w1, w2 = (x1 - x3) / (x6 - x1), (x2 - x3) / (x6 - x2) #???
    all(isfinite, (w1, w2)) && w1 != 0 && w2 != 0 && #???
        signbit(w1) == signbit(w2) || return nothing #???
    magnitude = sqrt(abs(w1)) * sqrt(abs(w2)) #???
    isfinite(magnitude) && magnitude > 0 || return nothing #???
    for scale in (magnitude, -magnitude) #???
        pole = scale == 1 ? Inf : (x3 - scale * x6) / (1 - scale) #???
        # !!! min(x3, x6) < pole < max(x3, x6) && continue
        #??? Полюс допустим между x₃ и x₆, если интервалы x₄…x₄′ и x₅…x₅′ его не пересекают.
        transform(s) = ((s - x3) - scale * (x6 - s)) / #???
                       ((s - x3) + scale * (x6 - s)) #???
        images = transform.(points) #???
        all(isfinite, images) || continue #???
        abs(images[1] + images[2]) <= #???
            1e-9 * max(1.0, abs(images[1]), abs(images[2])) || continue #???
        abs(images[3] + 1) <= 1e-9 && abs(images[6] - 1) <= 1e-9 || continue #???
        -1 < images[4] < images[5] < 1 || continue #???
        return (; transform, pole, gamma=images[1], images) #???
    end #???
    return nothing #???
end #???

function _gamma_add_primes(sample, s, t; N, M) #???
    x4prime = _gamma_finite_value(() -> x_4_prime(s, t; N=N, M=M)) #???
    x5prime = _gamma_finite_value(() -> x_5_prime(s, t; N=N, M=M)) #???
    # !!! all(isfinite, (x4prime, x5prime)) || return nothing
    #??? Дизъюнкция в условии допускает одну корректную prime-границу.
    return merge(sample, (; x4prime, x5prime)) #???
end #???

function _gamma_sample(s, t; N, M, primes=false) #???
    funcs = (x_1, x_2, x_3, x_4, x_5, x_6) #???
    points = ntuple(k -> _gamma_finite_value(() -> funcs[k](s, t; N=N, M=M)), 6) #???
    normalization = _gamma_normalization(points) #???
    normalization === nothing && return nothing #???
    sample = (; gamma=normalization.gamma, points, normalization) #???
    # !!! if !primes; return (; gamma=normalization.gamma, points, normalization); end
    return primes ? _gamma_add_primes(sample, s, t; N=N, M=M) : sample #???
end #???

#??? Корень γ(s,t)=γ₀ на ребре: секущая внутри скобки и бисекция при плохом шаге.
function _gamma_edge_sample(s1, t1, s2, t2, g1, g2, gamma0; N, M, maxiter=16) #???
    (g1 <= gamma0 < g2 || g2 <= gamma0 < g1) || return nothing #???
    low, high = 0.0, 1.0 #???
    f_low, f_high = g1 - gamma0, g2 - gamma0 #???
    λ = clamp(-f_low / (f_high - f_low), 0.0, 1.0) #???
    best = nothing #???
    for _ in 1:maxiter #???
        s, t = s1 + λ * (s2 - s1), t1 + λ * (t2 - t1) #???
        sample = _gamma_sample(s, t; N=N, M=M) #???
        sample === nothing && return nothing #???
        residual = sample.gamma - gamma0 #???
        # !!! best = (s, t, residual)
        best = (s, t, residual, sample) #???
        abs(residual) <= 1e-9 * max(1.0, abs(gamma0)) && break #???
        if (residual < 0) == (f_low < 0) #???
            low, f_low = λ, residual #???
        else #???
            high, f_high = λ, residual #???
        end #???
        proposed = low - f_low * (high - low) / (f_high - f_low) #???
        λ = isfinite(proposed) && low + 0.1 * (high - low) < proposed < #???
            high - 0.1 * (high - low) ? proposed : (low + high) / 2 #???
    end #???
    best === nothing && return nothing #???
    # !!! s, t, residual = best
    s, t, residual, sample = best #???
    abs(residual) <= 1e-5 * max(1.0, abs(gamma0)) || return nothing #???
    # !!! return _gamma_sample(s, t; N=N, M=M, primes=true)
    return _gamma_add_primes(sample, s, t; N=N, M=M) #??? Без повторного расчёта x₁…x₆.
end #???

function _gamma_clip(poly::Vector{_GammaPoint}, signed_distance) #???
    result = _GammaPoint[] #???
    for k in eachindex(poly) #???
        a, b = poly[k], poly[mod1(k + 1, length(poly))] #???
        fa, fb = signed_distance(a), signed_distance(b) #???
        inside_a, inside_b = fa > 0, fb > 0 #???
        if inside_a != inside_b #???
            weight = clamp(fa / (fa - fb), 0.0, 1.0) #???
            push!(result, (a[1] + weight * (b[1] - a[1]), #???
                           a[2] + weight * (b[2] - a[2]))) #???
        end #???
        inside_b && push!(result, b) #???
    end #???
    return result #???
end #???

function _gamma_contains(poly::Vector{_GammaPoint}, x, y) #???
    inside = false #???
    for k in eachindex(poly) #???
        a, b = poly[k], poly[mod1(k + 1, length(poly))] #???
        if (a[2] > y) != (b[2] > y) #???
            crossing = a[1] + (y - a[2]) * (b[1] - a[1]) / (b[2] - a[2]) #???
            x < crossing && (inside = !inside) #???
        end #???
    end #???
    return inside #???
end #???

function _gamma_polygon_area(poly) #???
    twice_area = 0.0 #???
    for k in eachindex(poly) #???
        a, b = poly[k], poly[mod1(k + 1, length(poly))] #???
        twice_area += a[1] * b[2] - b[1] * a[2] #???
    end #???
    return abs(twice_area) / 2 #???
end #???

function _gamma_paint!(mask, coordinates, triangle::Vector{_GammaPoint}; stats=nothing) #???
    all(p -> all(isfinite, p), triangle) || return 0 #???
    poly = triangle #???
    for distance in (p -> p[1] + 1, p -> 1 - p[1], #???
                     p -> p[2] + 1, p -> 1 - p[2], p -> p[2] - p[1]) #???
        poly = _gamma_clip(poly, distance) #???
        length(poly) >= 3 || return 0 #???
    end #???
    area = _gamma_polygon_area(poly) #???
    area > 0 || return 0 #???
    if stats !== nothing #???
        stats[1] += 1 #??? Число многоугольников с ненулевой площадью в квадрате.
        stats[2] += area #??? Сумма площадей до объединения.
    end #???
    xs, ys = extrema(p[1] for p in poly), extrema(p[2] for p in poly) #???
    i1, i2 = searchsortedfirst(coordinates, xs[1]), searchsortedlast(coordinates, xs[2]) #???
    j1, j2 = searchsortedfirst(coordinates, ys[1]), searchsortedlast(coordinates, ys[2]) #???
    painted = 0 #???
    for j in j1:j2, i in i1:i2 #???
        x, y = coordinates[i], coordinates[j] #???
        if !mask[i, j] && x < y && _gamma_contains(poly, x, y) #???
            mask[i, j] = true #???
            painted += 1 #???
        end #???
    end #???
    if painted == 0 #??? Слишком узкая область могла пройти между центрами пикселей.
        half_step = step(coordinates) / 2 #???
        ii1 = max(1, searchsortedfirst(coordinates, xs[1] - half_step)) #???
        ii2 = min(length(coordinates), searchsortedlast(coordinates, xs[2] + half_step)) #???
        jj1 = max(1, searchsortedfirst(coordinates, ys[1] - half_step)) #???
        jj2 = min(length(coordinates), searchsortedlast(coordinates, ys[2] + half_step)) #???
        for j in jj1:jj2, i in ii1:ii2 #???
            mask[i, j] && continue #???
            x, y = coordinates[i], coordinates[j] #???
            # !!! x - half_step < y + half_step || continue
            #??? Закрашиваемая точка сетки сама обязана удовлетворять -1<α₁<β₁<1.
            -1 < x < y < 1 || continue #???
            cell_part = poly #???
            for distance in (p -> p[1] - (x - half_step), #???
                             p -> (x + half_step) - p[1], #???
                             p -> p[2] - (y - half_step), #???
                             p -> (y + half_step) - p[2]) #???
                cell_part = _gamma_clip(cell_part, distance) #???
                length(cell_part) >= 3 || break #???
            end #???
            if length(cell_part) >= 3 && _gamma_polygon_area(cell_part) > 0 #???
                mask[i, j] = true #???
                painted += 1 #???
            end #???
        end #???
        painted > 0 && stats !== nothing && (stats[3] += 1) #???
    end #???
    return painted #???
end #???

#??? Полюс внутри маленького интервала по t требует разрыва, а не длинного ребра через ∞.
function _gamma_strip_crosses_pole(sample, family, t1, t2) #???
    p = sample.points #???
    a, b = family == 1 ? (sample.x4prime, p[4]) : (sample.x5prime, p[5]) #???
    v1, v2 = t1 * b + (1 - t1) * a, t2 * b + (1 - t2) * a #???
    return min(v1, v2) <= sample.normalization.pole <= max(v1, v2) #???
end #???

function _gamma_family_point(sample, family, t) #???
    p, T = sample.points, sample.normalization.transform #???
    if family == 1 #???
        return (T(t * p[4] + (1 - t) * sample.x4prime), T(p[5])) #???
    end #???
    return (T(p[4]), T(t * p[5] + (1 - t) * sample.x5prime)) #???
end #???

function _gamma_paint_segment!(mask, coordinates, a, b; nt=12, stats=nothing) #???
    painted = 0 #???
    for family in 1:2, k in 0:(nt - 1) #???
        if family == 1 #???
            isfinite(a.x4prime) && isfinite(b.x4prime) || continue #???
        else #???
            isfinite(a.x5prime) && isfinite(b.x5prime) || continue #???
        end #???
        t1, t2 = k / nt, (k + 1) / nt #???
        (_gamma_strip_crosses_pole(a, family, t1, t2) || #???
         _gamma_strip_crosses_pole(b, family, t1, t2)) && continue #???
        a1, b1 = _gamma_family_point(a, family, t1), #???
                 _gamma_family_point(b, family, t1) #???
        a2, b2 = _gamma_family_point(a, family, t2), #???
                 _gamma_family_point(b, family, t2) #???
        painted += _gamma_paint!(mask, coordinates, _GammaPoint[a1, b1, b2]; stats=stats) #???
        painted += _gamma_paint!(mask, coordinates, _GammaPoint[a1, b2, a2]; stats=stats) #???
    end #???
    return painted #???
end #???

""" #???
    gamma_value(x, y; N, M) #???

Вычисляет `γ=T(x₁)` по той же нормировке, что и `plot_gamma_slice`. #???
""" #???
function gamma_value(x, y; N::Int, M::Int) #???
    x > 0 && y > 0 || throw(ArgumentError("нужно x>0 и y>0")) #???
    1 <= M < N || throw(ArgumentError("нужно 1 ≤ M < N")) #???
    sample = _gamma_sample(x, y; N=N, M=M) #???
    if sample === nothing #???
        funcs = (x_1, x_2, x_3, x_4, x_5, x_6) #???
        anchors = ntuple(k -> _gamma_finite_value(() -> #???
            funcs[k](x, y; N=N, M=M)), 6) #???
        throw(ArgumentError("для x=$x, y=$y, N=$N, M=$M нормировка не определена; " * #???
                            "(x₁,x₂,x₃,x₄,x₅,x₆)=$anchors")) #???
    end #???
    return sample.gamma #???
end #???

""" #???
    find_gamma_reference(; N, M, x_bounds=(1e-4,50.0), y_bounds=(1e-8,50.0)) #???

Находит внутреннюю точку сетки с корректной нормировкой и пересечением #???
её уровня γ соседними узлами. Возвращает `(x, y, gamma)`. #???
""" #???
function find_gamma_reference(; N::Int, M::Int, #???
    # !!! x_bounds=(1e-3, 10.0), y_bounds=(1e-3, 10.0), nx=13, ny=13,
    x_bounds=(1e-4, 50.0), y_bounds=(1e-8, 50.0), nx=13, ny=13, #???
    preferred=(4.0, 0.026)) #???
    1 <= M < N || throw(ArgumentError("нужно 1 ≤ M < N")) #???
    nx >= 3 && ny >= 3 || throw(ArgumentError("nx и ny должны быть ≥3")) #???
    0 < x_bounds[1] < x_bounds[2] && 0 < y_bounds[1] < y_bounds[2] || #???
        throw(ArgumentError("диапазоны x и y должны быть положительными")) #???
    # !!! xs = exp.(range(log(x_bounds[1]), log(x_bounds[2]); length=nx))
    # !!! ys = exp.(range(log(y_bounds[1]), log(y_bounds[2]); length=ny))
    # !!! xs = _gamma_chebyshev_log_nodes(x_bounds, nx)
    # !!! ys = _gamma_chebyshev_log_nodes(y_bounds, ny)
    xs = exp.(range(log(x_bounds[1]), log(x_bounds[2]); length=nx)) #??? Равномерно по log(x).
    ys = exp.(range(log(y_bounds[1]), log(y_bounds[2]); length=ny)) #??? Равномерно по log(y).
    levels = fill(NaN, nx, ny) #???
    for j in 1:ny, i in 1:nx #???
        sample = _gamma_sample(xs[i], ys[j]; N=N, M=M) #???
        sample === nothing || (levels[i, j] = sample.gamma) #???
    end #???
    choices = Tuple{Float64, Int, Int}[] #???
    for j in 2:(ny - 1), i in 2:(nx - 1) #???
        value = levels[i, j] #???
        isfinite(value) || continue #???
        neighbors = (levels[i - 1, j], levels[i + 1, j], #???
                     levels[i, j - 1], levels[i, j + 1]) #???
        all(isfinite, neighbors) || continue #???
        any(v -> v < value, neighbors) && any(v -> v > value, neighbors) || continue #???
        distance = abs(log(xs[i] / preferred[1])) + #???
                   abs(log(ys[j] / preferred[2])) #???
        push!(choices, (distance, i, j)) #???
    end #???
    if isempty(choices) #???
        finite_count = count(isfinite, levels) #???
        throw(ArgumentError("не найдена внутренняя точка с пересекающимся уровнем γ; " * #???
                            "корректных узлов: $finite_count из $(nx * ny)")) #???
    end #???
    _, i, j = minimum(choices) #???
    return (; x=xs[i], y=ys[j], gamma=levels[i, j]) #???
end #???

#= !!! Отклонённая идея: искать треугольники по ступеням уже нарисованной маски.
#??? Ищем горизонтальные нижние границы чёрных ступеней.
#??? Их правые концы лежат около диагонали и задают вершины серых треугольников.
function _gamma_triangle_cuts(mask, coordinates; min_run_fraction=0.02, #???
                              tolerance_pixels=4) #???
    n = length(coordinates) #???
    lower = [findfirst(view(mask, i, :)) for i in 1:n] #???
    min_run = max(4, ceil(Int, min_run_fraction * n)) #???
    step_size = step(coordinates) #???
    candidates = Float64[] #???
    for level in 2:(n - 2) #???
        level_y = coordinates[level] #???
        run_start = 0 #???
        for i in 1:(n + 1) #???
            on_plateau = i <= n && lower[i] !== nothing && #???
                         abs(lower[i] - level) <= 1 #???
            if on_plateau && run_start == 0 #???
                run_start = i #???
            elseif !on_plateau && run_start != 0 #???
                last_i = i - 1 #???
                if last_i - run_start + 1 >= min_run && #???
                   abs(coordinates[last_i] - level_y) <= #???
                       tolerance_pixels * step_size #???
                    push!(candidates, coordinates[run_start], level_y) #???
                end #???
                run_start = 0 #???
            end #???
        end #???
    end #???
    isempty(candidates) && return Float64[] #???
    sort!(candidates) #???
    cuts = Float64[-1.0] #???
    for value in candidates #???
        value > cuts[end] + tolerance_pixels * step_size && #???
            value < 1 - tolerance_pixels * step_size && #???
            push!(cuts, value) #???
    end #???
    push!(cuts, 1.0) #???
    return length(cuts) >= 3 ? cuts : Float64[] #???
end #???

function _gamma_manual_triangle_cuts(values) #???
    cuts = sort!(Float64.(collect(values))) #???
    all(isfinite, cuts) && all(v -> -1 <= v <= 1, cuts) || #???
        throw(ArgumentError("вершины треугольников должны лежать в [-1,1]")) #???
    isempty(cuts) || first(cuts) > -1 && pushfirst!(cuts, -1.0) #???
    isempty(cuts) || last(cuts) < 1 && push!(cuts, 1.0) #???
    unique!(cuts) #???
    length(cuts) >= 2 || throw(ArgumentError("нужно хотя бы две разные вершины")) #???
    return cuts #???
end #???
=#

function _gamma_triangle_vertices(sample) #???
    isfinite(sample.x4prime) && isfinite(sample.x5prime) || return nothing #???
    p, T = sample.points, sample.normalization.transform #???
    vertices = _GammaPoint[(T(p[4]), T(p[5])), #???
                           (T(sample.x4prime), T(p[5])), #???
                           (T(p[4]), T(sample.x5prime))] #???
    all(v -> all(isfinite, v), vertices) || return nothing #???
    return vertices #???
end #???

#??? Треугольник с заданными пользователем вершинами перемещается вдоль
#??? короткого участка линии γ=γ₀. Берём сечения самого треугольника.
function _gamma_paint_triangle_sweep!(mask, coordinates, a, b; #???
                                      steps=7, stats=nothing) #???
    av, bv = _gamma_triangle_vertices(a), _gamma_triangle_vertices(b) #???
    if av === nothing || bv === nothing #???
        return 0 #???
    end #???
    #= !!! Прежнее заполнение полос между рёбрами треугольников.
    painted = _gamma_paint!(mask, coordinates, av; stats=stats) +
              _gamma_paint!(mask, coordinates, bv; stats=stats)
    for (i, j) in ((1, 2), (2, 3), (3, 1))
        painted += _gamma_paint!(mask, coordinates,
                                 _GammaPoint[av[i], av[j], bv[j]]; stats=stats)
        painted += _gamma_paint!(mask, coordinates,
                                 _GammaPoint[av[i], bv[j], bv[i]]; stats=stats)
    end
    =#
    painted = 0 #???
    for λ in range(0.0, 1.0; length=steps) #???
        triangle = _GammaPoint[ntuple(k -> (1 - λ) * av[i][k] + #???
                                         λ * bv[i][k], 2) for i in 1:3] #???
        painted += _gamma_paint!(mask, coordinates, triangle; stats=stats) #???
    end #???
    return painted #???
end #???

#??? Концы горизонтальных ступеней разбивают всю диагональ от -1 до 1.
#??? Каждая соседняя пара a,b задаёт треугольник (a,a), (b,b), (a,b).
function _gamma_diagonal_intervals(mask, coordinates; #???
                                   min_run_fraction=0.02, tolerance_pixels=4) #???
    n = length(coordinates) #???
    lower = [findfirst(view(mask, i, :)) for i in 1:n] #???
    min_run = max(4, ceil(Int, min_run_fraction * n)) #???
    tolerance = tolerance_pixels * step(coordinates) #???
    candidates = Tuple{Float64, Float64}[] #???
    for level in 2:(n - 1) #???
        run_start = 0 #???
        for i in 1:(n + 1) #???
            on_level = i <= n && lower[i] !== nothing && #???
                       abs(lower[i] - level) <= 1 #???
            if on_level && run_start == 0 #???
                run_start = i #???
            elseif !on_level && run_start != 0 #???
                last_i = i - 1 #???
                a, b = coordinates[run_start], coordinates[level] #???
                if last_i - run_start + 1 >= min_run && #???
                   b > a + tolerance && #???
                   abs(coordinates[last_i] - b) <= tolerance #???
                    push!(candidates, (a, b)) #???
                end #???
                run_start = 0 #???
            end #???
        end #???
    end #???
    sort!(candidates) #???
    intervals = Tuple{Float64, Float64}[] #???
    for (a, b) in candidates #???
        any(old -> abs(old[1] - a) <= tolerance && #???
                   abs(old[2] - b) <= tolerance, intervals) && continue #???
        push!(intervals, (a, b)) #???
    end #???
    isempty(intervals) && return intervals #??? Без надёжных ступеней разбиение неизвестно.
    cuts = Float64[-1.0] #???
    for (a, b) in intervals, value in (a, b) #???
        value > cuts[end] + tolerance && value < 1 - tolerance && #???
            push!(cuts, value) #???
    end #???
    push!(cuts, 1.0) #???
    length(cuts) >= 3 || return Tuple{Float64, Float64}[] #???
    return [(cuts[k], cuts[k + 1]) for k in 1:(length(cuts) - 1)] #???
end #???

function _gamma_reflect_antidiagonal!(classes) #???
    n = size(classes, 1) #???
    size(classes, 2) == n || throw(ArgumentError("для отражения нужна квадратная сетка")) #???
    for j in 2:n, i in 1:(j - 1) #???
        if i + j < n + 1 #??? Точка ниже β₁=-α₁; источник выше неё.
            classes[i, j] = classes[n + 1 - j, n + 1 - i] #???
        end #???
    end #???
    return classes #???
end #???

""" #???
    plot_gamma_slice(gamma0; N, M=nothing, x_bounds=(1e-4,50.0), #???
                     y_bounds=(1e-8,50.0), nx=50, ny=50, pixels=800) #???

При фиксированном `γ=T(x₁)` численно закрашивает в квадрате `[-1,1]²` #???
точки из двух заданных семейств с `α₁<β₁`. `M=nothing` объединяет все #???
`M=1:floor(N/2)`. #???
Концы нижних горизонтальных ступеней разбивают диагональ на отрезки; #???
каждый отрезок от `a` до `b` задаёт серый треугольник с вершинами #???
`(a,a)`, `(b,b)`, `(a,b)`. #???
Все треугольники имеют один оттенок; остальная часть над диагональю — другой. #???
Часть с `α₁+β₁>0` отражается относительно прямой `β₁=-α₁`. #???
Узлы `(nx,ny)` равномерны по `log(x),log(y)` и включают концы диапазонов. #???
Точность определяется сеткой `(nx,ny)`, `nt` и `pixels`. #???
""" #???
# !!! Треугольники с вершинами `(T(x₄),T(x₅))`, `(T(x₄′),T(x₅))`,
# !!! `(T(x₄),T(x₅′))` перекрашивали почти всю чёрную область.
function plot_gamma_slice(gamma0; N::Int, M=nothing, #???
    # !!! x_bounds=(1e-3, 10.0), y_bounds=(1e-3, 10.0),
    x_bounds=(1e-4, 50.0), y_bounds=(1e-8, 50.0), #???
    nx::Int=50, ny::Int=50, pixels::Int=800, nt::Int=12, #???
    # !!! triangle_cuts=nothing,
    # !!! triangle_shades=(0.76, 0.59, 0.86, 0.66, 0.94, 0.70),
    triangle_shade::Real=0.65, remaining_shade::Real=0.90, #???
    dark_shade::Real=0.16, #??? Бывший чёрный цвет стал тёмно-серым.
    graph_name::AbstractString="gamma_slice") #???
    isfinite(gamma0) || throw(ArgumentError("γ должно быть конечным")) #???
    N >= 2 || throw(ArgumentError("N должно быть не меньше 2")) #???
    0 < x_bounds[1] < x_bounds[2] && 0 < y_bounds[1] < y_bounds[2] || #???
        throw(ArgumentError("диапазоны x и y должны иметь положительные возрастающие концы")) #???
    nx >= 2 && ny >= 2 && pixels >= 2 && nt >= 1 || #???
        throw(ArgumentError("nx, ny, pixels должны быть ≥2, а nt ≥1")) #???
    # !!! !isempty(triangle_shades) && all(v -> isfinite(v) && 0 <= v <= 1, triangle_shades) ||
    # !!!     throw(ArgumentError("triangle_shades должен содержать оттенки от 0 до 1"))
    # !!! all(v -> isfinite(v) && 0 <= v <= 1, (triangle_shade, remaining_shade)) &&
    all(v -> isfinite(v) && 0 <= v <= 1, #???
        (triangle_shade, remaining_shade, dark_shade)) && #???
        triangle_shade != remaining_shade || #???
        throw(ArgumentError("оттенки должны лежать между 0 и 1, а два светлых оттенка различаться")) #???
    isempty(strip(graph_name)) && throw(ArgumentError("graph_name не может быть пустым")) #???
    # !!! ms = isnothing(M) ? collect(1:(N - 1)) :
    ms = isnothing(M) ? collect(1:fld(N, 2)) : #???
         (M isa Integer ? [M] : collect(M)) #???
    isempty(ms) && throw(ArgumentError("список M пуст")) #???
    # !!! all(m -> m isa Integer && 1 <= m < N, ms) ||
    all(m -> m isa Integer && 1 <= m <= fld(N, 2), ms) || #???
        throw(ArgumentError("каждое M должно лежать в 1:floor(N/2)")) #???
    ms = sort!(unique(Int.(ms))) #???
    coordinates = range(-1 + 1 / pixels, 1 - 1 / pixels; length=pixels) #???
    mask = falses(pixels, pixels) #???
    # !!! triangle_masks = BitMatrix[]
    gamma0 = Float64(gamma0) #???
    # !!! sx = _gamma_chebyshev_log_nodes(x_bounds, nx)
    # !!! sy = _gamma_chebyshev_log_nodes(y_bounds, ny)
    sx = exp.(range(log(x_bounds[1]), log(x_bounds[2]); length=nx)) #??? Равномерно по log(x).
    sy = exp.(range(log(y_bounds[1]), log(y_bounds[2]); length=ny)) #??? Равномерно по log(y).
    for m in ms #???
        started = time_ns() #???
        # !!! sx = exp.(range(log(x_bounds[1]), log(x_bounds[2]); length=nx))
        # !!! sy = exp.(range(log(y_bounds[1]), log(y_bounds[2]); length=ny))
        levels = fill(NaN, nx, ny) #???
        # !!! progress_step = max(1, cld(ny, 5))
        progress_step = max(1, cld(ny, 20)) #??? Чаще показываем ход расчёта при расширенном диапазоне.
        first_row_step = max(1, cld(nx, 10)) #??? Первый ряд может быть самым дорогим при малом y.
        for j in 1:ny #???
            for i in 1:nx #???
                if j == 1 && (i == 1 || i % first_row_step == 0) #???
                    @info "Первый ряд сетки γ" M=m узел=i всего=nx x=sx[i] y=sy[j] секунд=round((time_ns() - started) / 1e9; digits=1) #???
                end #???
                sample = _gamma_sample(sx[i], sy[j]; N=N, M=m) #???
                sample === nothing || (levels[i, j] = sample.gamma) #???
            end #???
            # !!! if j % progress_step == 0 || j == ny
            if j == 1 || j % progress_step == 0 || j == ny #???
                # !!! @info "Вычисление γ" M=m строк=j всего=ny
                @info "Вычисление γ" M=m строк=j всего=ny секунд=round((time_ns() - started) / 1e9; digits=1) #???
            end #???
        end #???
        grid_seconds = round((time_ns() - started) / 1e9; digits=1) #???
        finite_levels = filter(isfinite, vec(levels)) #???
        level_range = isempty(finite_levels) ? nothing : extrema(finite_levels) #???
        # !!! @info "Сетка γ" M=m конечных=length(finite_levels) диапазон=level_range цель=gamma0
        @info "Сетка γ" M=m конечных=length(finite_levels) диапазон=level_range цель=gamma0 секунд=grid_seconds кеш_тета=length(THETA_MATRIX_CACHE) кеш_преобразований=length(MODULAR_TRANSFORM_CACHE) кеш_характеристик=length(MODULAR_CHAR_CACHE) #???
        contour_started = time_ns() #???
        edge_cache = Dict{Tuple{Tuple{Int, Int}, Tuple{Int, Int}}, Any}() #???
        edge_time_ns = Ref(0) #??? Время вычисления пересечений уровня γ с рёбрами.
        edge_count = Ref(0) #??? Число действительно вычисленных рёбер без попаданий в кеш.
        function crossing(p, q) #???
            key = p < q ? (p, q) : (q, p) #???
            return get!(edge_cache, key) do #???
                g1, g2 = levels[p...], levels[q...] #???
                # !!! _gamma_edge_sample(sx[p[1]], sy[p[2]], sx[q[1]], sy[q[2]],
                # !!!                    g1, g2, gamma0; N=N, M=m)
                edge_started = time_ns() #???
                result = _gamma_edge_sample(sx[p[1]], sy[p[2]], sx[q[1]], sy[q[2]], #???
                                            g1, g2, gamma0; N=N, M=m) #???
                edge_time_ns[] += time_ns() - edge_started #???
                edge_count[] += 1 #???
                result #???
            end #???
        end #???
        segments, painted = 0, 0 #???
        paint_time_ns = 0 #??? Время построения и закраски полос на сетке пикселей.
        # !!! triangle_mask = falses(pixels, pixels)
        # !!! triangle_painted = 0
        first_available, second_available = 0, 0 #???
        paint_stats = [0.0, 0.0, 0.0] #??? Многоугольники, площадь, субпиксельные попадания.
        # !!! triangle_stats = [0.0, 0.0, 0.0]
        for j in 1:(ny - 1), i in 1:(nx - 1) #???
            for corners in (((i, j), (i + 1, j), (i + 1, j + 1)), #???
                            ((i, j), (i + 1, j + 1), (i, j + 1))) #???
                all(p -> isfinite(levels[p...]), corners) || continue #???
                gmin, gmax = extrema(levels[p...] for p in corners) #???
                gmin <= gamma0 < gmax || continue #???
                cuts = Any[] #???
                for (a, b) in ((1, 2), (2, 3), (3, 1)) #???
                    sample = crossing(corners[a], corners[b]) #???
                    sample === nothing || push!(cuts, sample) #???
                end #???
                length(cuts) == 2 || continue #???
                segments += 1 #???
                first_available += all(isfinite(cut.x4prime) for cut in cuts) #???
                second_available += all(isfinite(cut.x5prime) for cut in cuts) #???
                # !!! painted += _gamma_paint_segment!(mask, coordinates, cuts[1], cuts[2];
                # !!!                                   nt=nt, stats=paint_stats)
                paint_started = time_ns() #???
                painted += _gamma_paint_segment!(mask, coordinates, cuts[1], cuts[2]; #???
                                                  nt=nt, stats=paint_stats) #???
                paint_time_ns += time_ns() - paint_started #???
                # !!! triangle_painted += _gamma_paint_triangle_sweep!(triangle_mask,
                # !!!     coordinates, cuts[1], cuts[2]; stats=triangle_stats)
            end #???
        end #???
        # !!! push!(triangle_masks, triangle_mask)
        elapsed = round((time_ns() - started) / 1e9; digits=1) #???
        contour_seconds = round((time_ns() - contour_started) / 1e9; digits=1) #???
        # !!! @info "Готово для M" M=m сегментов=segments первая_граница=first_available вторая_граница=second_available новых_пикселей=painted многоугольников=Int(paint_stats[1]) площадь_до_объединения=paint_stats[2] субпиксельных=Int(paint_stats[3]) секунд=elapsed
        # !!! @info "Готово для M" M=m сегментов=segments первая_граница=first_available вторая_граница=second_available новых_пикселей=painted многоугольников=Int(paint_stats[1]) площадь_до_объединения=paint_stats[2] субпиксельных=Int(paint_stats[3]) сетка_секунд=grid_seconds пересечения_и_заливка_секунд=contour_seconds секунд=elapsed
        @info "Готово для M" M=m сегментов=segments первая_граница=first_available вторая_граница=second_available новых_пикселей=painted многоугольников=Int(paint_stats[1]) площадь_до_объединения=paint_stats[2] субпиксельных=Int(paint_stats[3]) сетка_секунд=grid_seconds пересечений=edge_count[] пересечения_секунд=round(edge_time_ns[] / 1e9; digits=1) заливка_секунд=round(paint_time_ns / 1e9; digits=1) пересечения_и_заливка_секунд=contour_seconds секунд=elapsed #???
        # !!! @info "Треугольники для M" M=m пикселей=count(triangle_mask) новых_пикселей=triangle_painted площадь_до_объединения=triangle_stats[2]
    end #???
    any(mask) || @warn "Чёрная область пуста на выбранной сетке; проверьте диапазоны γ в сообщениях выше" #???
    diagonal_intervals = _gamma_diagonal_intervals(mask, coordinates) #???
    isempty(diagonal_intervals) && @warn "Не найдены ступени для построения диагональных треугольников; попробуйте увеличить pixels" #???
    @info "Треугольники на диагонали" интервалы=diagonal_intervals #???
    fig = Makie.Figure(size=(1100, 1050)) #???
    # !!! ax = Makie.Axis(fig[1, 1]; xlabel="α₁", ylabel="β₁", aspect=Makie.DataAspect())
    ax = Makie.Axis(fig[1, 1]; xlabel="", ylabel="", #???
                    aspect=Makie.DataAspect()) #???
    #= !!! Несколько слоёв heatmap скрывали треугольники или перекрашивали чёрное непредсказуемо.
    for (k, triangle_mask) in enumerate(triangle_masks)
        visible_triangle = triangle_mask .& .!mask
        any(visible_triangle) || continue
        shade = triangle_shades[mod1(k, length(triangle_shades))]
        Makie.heatmap!(ax, coordinates, coordinates, Float32.(visible_triangle);
                       colormap=[Makie.RGBAf(1, 1, 1, 0),
                                 Makie.RGBAf(shade, shade, shade, 1)],
                       colorrange=(0, 1))
    end
    Makie.heatmap!(ax, coordinates, coordinates, Float32.(mask);
                   colormap=[Makie.RGBAf(1, 1, 1, 0), Makie.RGBAf(0, 0, 0, 1)],
                   colorrange=(0, 1))
    =#
    #= !!! Категориальная заливка серым перекрывала прежнюю чёрную фигуру.
    classes = zeros(Float32, pixels, pixels)
    classes[mask] .= 1
    colors = Makie.RGBAf[Makie.RGBAf(1, 1, 1, 1), Makie.RGBAf(0, 0, 0, 1)]
    for (k, triangle_mask) in enumerate(triangle_masks)
        shade = triangle_shades[mod1(k, length(triangle_shades))]
        push!(colors, Makie.RGBAf(shade, shade, shade, 1))
        classes[triangle_mask] .= Float32(k + 1)
    end
    Makie.heatmap!(ax, coordinates, coordinates, classes;
                   colormap=colors, colorrange=(0, length(colors) - 1))
    =#
    # !!! Makie.heatmap!(ax, coordinates, coordinates, Float32.(mask);
    # !!!                colormap=[Makie.RGBAf(1, 1, 1, 1), Makie.RGBAf(0, 0, 0, 1)],
    # !!!                colorrange=(0, 1))
    # !!! classes = zeros(Float32, pixels, pixels)
    # !!! classes[mask] .= 1
    # !!! colors = Makie.RGBAf[Makie.RGBAf(1, 1, 1, 1), Makie.RGBAf(0, 0, 0, 1)]
    classes = zeros(Float32, pixels, pixels) #??? 0 — белое под диагональю.
    for j in 2:pixels #???
        classes[1:(j - 1), j] .= 1 #??? 1 — светло-серое над диагональю.
    end #???
    # !!! classes[mask] .= 2
    classes[mask] .= 2 #??? 2 — тёмно-серая область.
    colors = Makie.RGBAf[Makie.RGBAf(1, 1, 1, 1), #???
                         Makie.RGBAf(remaining_shade, remaining_shade, remaining_shade, 1), #???
                         # !!! Makie.RGBAf(0, 0, 0, 1),
                         Makie.RGBAf(dark_shade, dark_shade, dark_shade, 1), #???
                         Makie.RGBAf(triangle_shade, triangle_shade, triangle_shade, 1)] #???
    # !!! for (k, (a, b)) in enumerate(diagonal_intervals)
    for (a, b) in diagonal_intervals #???
        triangle_mask = falses(pixels, pixels) #???
        vertices = _GammaPoint[(a, a), (b, b), (a, b)] #???
        _gamma_paint!(triangle_mask, coordinates, vertices) #???
        # !!! shade = triangle_shades[mod1(k, length(triangle_shades))]
        # !!! push!(colors, Makie.RGBAf(shade, shade, shade, 1))
        # !!! classes[triangle_mask] .= Float32(k + 1)
        visible_triangle = triangle_mask .& .!mask #??? Сохраняем прежнюю чёрную фигуру.
        # !!! classes[visible_triangle] .= Float32(k + 1)
        classes[visible_triangle] .= 3 #??? Один оттенок для всех треугольников.
        @info "Диагональный треугольник" вершины=vertices пикселей=count(visible_triangle) пересечение_с_чёрным=count(triangle_mask .& mask) #???
    end #???
    #??? Отражение (α₁,β₁) ↦ (-β₁,-α₁). Отражаем готовую раскраску,
    #??? чтобы границы тёмной области и треугольники совпали попиксельно.
    _gamma_reflect_antidiagonal!(classes) #???
    Makie.heatmap!(ax, coordinates, coordinates, classes; #???
                   colormap=colors, colorrange=(0, length(colors) - 1)) #???
    # !!! Если рисовать контур исходной чёрной маски поверх серой,
    # !!! внутри перекрашенных треугольников остаются лишние чёрные линии.
    # !!! if any(mask)
    # !!!     Makie.contour!(ax, coordinates, coordinates, Float32.(mask);
    # !!!                    levels=[0.5], color=:black, linewidth=1)
    # !!! end
    Makie.lines!(ax, [-1.0, 1.0], [-1.0, 1.0]; color=:gray, #???
                 linestyle=:dash, linewidth=1) #???
    Makie.xlims!(ax, -1, 1) #???
    Makie.ylims!(ax, -1, 1) #???
    file = endswith(lowercase(graph_name), ".png") ? graph_name : graph_name * ".png" #???
    path = joinpath(@__DIR__, file) #???
    Makie.save(path, fig) #???
    println("График сохранён: ", path) #???
    # !!! display(fig)
    if lowercase(get(ENV, "SHTIFFEL_PLOT_BACKEND", "gl")) != "cairo" #???
        display(fig) #??? В оконном режиме показываем рисунок как прежде.
    end #??? При CairoMakie достаточно уже сохранённого PNG.
    return fig #???
end #???
