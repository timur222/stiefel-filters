# Numerical inversion of the two integer height equations, including all
# boundary faces where one or both hidden prevertices coincide with a given end.

function _period_error(g,n,m1,m2)
    # Extra bits for subtraction in tile units, using precisely the stored
    # quadrature values. Accuracy of those values is checked by order doubling.
    a,b,c = BigFloat.(g.H)
    total = a+b+c
    Float64[n*a/total-m1,n*(a+b)/total-m2]
end

function _mode_candidates(E,n,order)
    coords(g) = n .* cumsum(g.H)[1:2] ./ sum(g.H)
    coarse, q = coords(geometry(E,order)),coords(geometry(E,2order))
    norm(q-coarse,Inf) <= 0.5 || throw(GeometryError("unstable mode window; raise quadrature order"))
    # The analytic flux estimate is ±2 tiles; one extra tile absorbs rounding.
    ranges = (max(1,ceil(Int,q[1]-3)):min(n-2,floor(Int,q[1]+3)),
              max(2,ceil(Int,q[2]-3)):min(n-1,floor(Int,q[2]+3)))
    pairs = [(a,b) for a in ranges[1] for b in ranges[2] if a<b]
    sort!(pairs, by=t->(hypot(t[1]-q[1],t[2]-q[2]),t))
end

_mouth(side,j) = 2j + (side == :below ? 2 : 3)
_trim_wall(side,j) = j + (side == :below ? 1 : 2)

function _geometry_on_face(E,free,v,sides,order)
    e = copy(E)
    for (k,j) in enumerate(free)
        i = 2j+2
        s = 1/(1+exp(-v[k]))
        e[_mouth(sides[j],j)] = E[i]+(E[i+1]-E[i])*s
    end
    geometry(e,order)
end

_recoverable(err) = err isa GeometryError || err isa SingularException ||
                    err isa DomainError || err isa PosDefException || err isa ArgumentError

function _solve_face(E,free,levels,sides,n,order)
    m1,m2 = levels
    function at(v)
        g = _geometry_on_face(E,free,v,sides,order)
        g,_period_error(g,n,m1,m2)
    end
    v = zeros(length(free))
    g,r = at(v)
    tol = min(n*2e-11,1e-8)
    # A mouth very close to a given endpoint can move the inferred trim
    # substantially. Refine to the attainable tile precision before accepting.
    refine_tol = min(tol,max(2e-14,8eps(Float64)*n))
    for step in 0:100
        norm(r,Inf) <= refine_tol && return g
        isempty(free) && break
        J = Matrix{Float64}(undef,2,length(v))
        for k in eachindex(v)
            d = zeros(length(v))
            d[k] = 0.002
            J[:,k] = (last(at(v+d))-last(at(v-d)))/(2d[k])
        end
        # SVD also handles one free mouth: J is then a two-by-one matrix.
        direction = -(pinv(J)*r)
        all(isfinite,direction) || break
        norm(direction,Inf) > 6 && (direction .*= 6/norm(direction,Inf))
        accepted = false
        for j in 0:15
            trial = clamp.(v + 0.5^j*direction,-28.,28.)
            trial == v && continue
            try
                newg,newr = at(trial)
                if norm(newr) < norm(r)
                    v,g,r = trial,newg,newr
                    accepted = true
                    break
                end
            catch err
                _recoverable(err) || rethrow()
            end
        end
        accepted || break
    end
    norm(r,Inf) <= tol && return g
    throw(GeometryError("period solve failed on face of dimension $(length(free)) (tile residual $(norm(r,Inf)))"))
end

function _check_heights(g,m,T)
    a,b,c = BigFloat.(g.H)
    Tbig = BigFloat(T)
    absmax = maximum(abs.([a,a+b,a+b+c,BigFloat(g.left)] ./ Tbig .- [m.m1,m.m2,m.n,m.n]))
    all(x->-1<x<1,g.slit) || throw(GeometryError("slit tip outside the open rectangle"))
    Float64(absmax)
end

function _recover_end(g,band,height)
    a,b = g.e[2band-1:2band]
    for _ in 1:60
        middle = (a+b)/2
        middle in (a,b) && break
        if _wall_height(g,middle,band) < height
            a = middle
        else
            b = middle
        end
    end
    (a+b)/2
end

function _validate_candidate(g,E,m1,m2,n,sides,free; T=sum(g.H)/n, saved_trims=nothing)
    fine = geometry(g.e,2g.order)
    degree = (;m1,m2,n)
    height_error = max(_check_heights(g,degree,T),_check_heights(fine,degree,T))
    height_error <= 1e-7 || throw(GeometryError("wall height error $(height_error) tiles at doubled order"))
    levels = (m1,m2)
    trims = zeros(2)
    endpoint_error = 0.0
    trim_error = 0.0
    for j in 1:2
        i,band = _mouth(sides[j],j),_trim_wall(sides[j],j)
        x = E[i]
        sign = sides[j] == :below ? -1.0 : 1.0
        measured = sign * Float64(BigFloat(_wall_height(g,x,band))/T-levels[j])
        delta = isnothing(saved_trims) ? (j in free ? measured : 0.0) : saved_trims[j]
        -1e-10 <= delta <= 1+1e-10 || throw(GeometryError("trim outside [0,1] at slit $j"))
        trims[j] = clamp(delta,0.,1.)
        target = BigFloat(levels[j]) + sign*BigFloat(trims[j])
        error = abs(BigFloat(_wall_height(fine,x,band))/T-target)
        trim_error = max(trim_error,Float64(error))
        if trims[j] != 0
            rebuilt = _recover_end(fine,band,Float64(target*T))
            endpoint_error = max(endpoint_error,abs(rebuilt-x))
        else
            endpoint_error = max(endpoint_error,abs(g.e[i]-x))
        end
    end
    trim_error <= 1e-7 || throw(GeometryError("trim error $(trim_error) tiles at doubled order"))
    endpoint_error <= max(1e-9*(last(E)-first(E)),8eps(maximum(abs,E))) ||
        throw(GeometryError("endpoint reconstruction discrepancy $(endpoint_error)"))
    lengths = Float64[m1,m2-m1,n-m2]
    for j in 1:2
        lengths[_trim_wall(sides[j],j)-1] -= trims[j]
    end
    minimum(lengths)>0 || throw(GeometryError("one of the trimmed passbands collapsed"))
    (T=T,trims=Tuple(trims),max_height_error_tiles=height_error,
     max_trim_error_tiles=trim_error,endpoint_error=endpoint_error)
end

"""
    synthesize(endpoints, n; levels=nothing, sides=nothing, order=128)

Construct the first numerically checked four-band mode. Endpoints are strictly
ordered, with the first pair delimiting the stopband. Optionally fix two
integer slit levels and/or two trim sides (:above or :below).
"""
function synthesize(endpoints,n::Integer;levels=nothing,sides=nothing,order::Int=128)
    n isa Bool && throw(ArgumentError("expected integer degree"))
    3 <= n < 2^53 || throw(ArgumentError("degree must be between 3 and 2^53-1"))
    2 <= order || throw(ArgumentError("quadrature order must be at least two"))
    E = Float64.(collect(endpoints))
    _check_endpoints(E)
    if levels === nothing
        pairs = _mode_candidates(E,Int(n),order)
    else
        length(levels)==2 && all(m->m isa Integer && !(m isa Bool),levels) &&
            0<levels[1]<levels[2]<n || throw(ArgumentError("require integer 0<m1<m2<n"))
        pairs = [(Int(levels[1]),Int(levels[2]))]
    end
    if sides === nothing
        side_pairs = [(a,b) for a in (:above,:below) for b in (:above,:below)]
    else
        length(sides)==2 && all(s->s in (:above,:below,"above","below"),sides) ||
            throw(ArgumentError("expected two trim sides, above or below"))
        side_pairs = [(Symbol(sides[1]),Symbol(sides[2]))]
    end
    attempts = 0
    last_failure = "no candidate levels"
    for (m1,m2) in pairs, choices in side_pairs
        attempts += 1
        for free in ((),(1,),(2,),(1,2))
            try
                g = _solve_face(E,free,(m1,m2),choices,n,order)
                checks = _validate_candidate(g,E,m1,m2,n,choices,free)
                diagnostics = (attempts=attempts,candidate_count=length(pairs)*length(side_pairs),
                    free_prevertex_count=length(free),max_height_error_tiles=checks.max_height_error_tiles,
                    max_trim_error_tiles=checks.max_trim_error_tiles,endpoint_error=checks.endpoint_error)
                return FilterModel(g,E,Int(n),m1,m2,checks.T,checks.trims,choices,diagnostics)
            catch err
                _recoverable(err) || rethrow()
                last_failure = sprint(showerror,err)
            end
        end
    end
    throw(GeometryError("no mode passed numerical checks ($attempts tried): $last_failure"))
end

"""Restore a saved compact model and recheck periods, trims and working endpoints."""
function restore_model(;e,endpoints,n,m1,m2,sides,trims,T,order=128)
    all(m->m isa Integer && !(m isa Bool),(n,m1,m2)) && 0<m1<m2<n<2^53 ||
        throw(ArgumentError("inconsistent saved integer mode"))
    length(sides)==2 && all(s->s in (:above,:below,"above","below"),sides) ||
        throw(ArgumentError("invalid saved trim sides"))
    length(trims)==2 && all(d->isfinite(d) && 0<=d<=1,trims) ||
        throw(ArgumentError("invalid saved trims"))
    isfinite(T) && T>0 || throw(ArgumentError("invalid saved period"))
    E = Float64.(collect(endpoints))
    _check_endpoints(E)
    g = geometry(e,order)
    choices = (Symbol(sides[1]),Symbol(sides[2]))
    indices = [_mouth(choices[j],j) for j in 1:2]
    all(i -> i in indices || E[i]==g.e[i], 1:8) ||
        throw(ArgumentError("unchanged working endpoints disagree with saved prevertices"))
    for b in 1:4
        g.e[2b-1]<=E[2b-1]<E[2b]<=g.e[2b] ||
            throw(ArgumentError("saved working band lies outside its maximal band"))
    end
    free = [j for j in 1:2 if E[indices[j]]!=g.e[indices[j]]]
    checks = _validate_candidate(g,E,m1,m2,n,choices,free; T=Float64(T),saved_trims=trims)
    diagnostics = (attempts=0,candidate_count=0,free_prevertex_count=length(free),
        max_height_error_tiles=checks.max_height_error_tiles,max_trim_error_tiles=checks.max_trim_error_tiles,
        endpoint_error=checks.endpoint_error)
    FilterModel(g,E,n,m1,m2,checks.T,checks.trims,choices,diagnostics)
end
