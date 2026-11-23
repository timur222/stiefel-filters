# The three projective normalization conditions e₁=-1, e₂=0, e₈=1 leave
# five free prevertices. Softmax coordinates keep the six gaps positive.
function _polygon_points(v)
    logits = [v;0.]
    weights = exp.(logits .- maximum(logits))
    cumulative = cumsum(weights)/sum(weights)
    [-1.;0.;cumulative[1:end-1];1.]
end

function _polygon_seed(initial)
    isnothing(initial) && return zeros(5)
    e = Float64.(collect(initial))
    _check_endpoints(e)
    a,b,c = e[1],e[2],e[8]
    gauge(x) = (c-a)*(x-b)/(2*(c-b)*(x-a)-(c-a)*(x-b))
    gaps = diff(gauge.(e[2:8]))
    all(>(0),gaps) || throw(ArgumentError("invalid polygon initialization"))
    log.(gaps[1:5]/gaps[6])
end

"""
    solve_polygon(n,m1,m2,T,h1,h2; trims=(0.,0.), sides=(:above,:above), initial=nothing, order=128)

Prescribe the two horizontal slit-tip abscissas and elliptic period. Solve
five equations in five prevertices in the projective gauge (-1,0,1), and
optionally trim the resulting working bands by the requested tile fractions.
"""
function solve_polygon(n::Integer,m1::Integer,m2::Integer,T::Real,h1::Real,h2::Real;
                       trims=(0.,0.),sides=(:above,:above),initial=nothing,order::Int=128)
    0<m1<m2<n<2^53 && all(x->!(x isa Bool),(n,m1,m2)) ||
        throw(ArgumentError("require 0<m1<m2<n<2^53"))
    isfinite(T) && T>0 && -1<h1<1 && -1<h2<1 ||
        throw(ArgumentError("expected a positive period and slit abscissas in (-1,1)"))
    length(trims)==2 && all(d->isfinite(d) && 0<=d<=1,trims) ||
        throw(ArgumentError("two trim lengths between zero and one required"))
    length(sides)==2 && all(s->s in (:above,:below,"above","below"),sides) ||
        throw(ArgumentError("two trim sides required"))
    choices = (Symbol(sides[1]),Symbol(sides[2]))
    counts = [m1,m2-m1,n-m2]
    targets = T .* counts
    retained = Float64.(counts)
    for j in 1:2
        retained[_trim_wall(choices[j],j)-1] -= trims[j]
    end
    minimum(retained)>0 || throw(ArgumentError("a retained passband would collapse"))
    function at(v)
        g = geometry(_polygon_points(v),order)
        g,[log.(g.H ./ targets);g.slit.-[h1,h2]]
    end
    v = _polygon_seed(initial)
    g,r = at(v)
    iterations = 0
    while norm(r)>2e-12 && iterations<80
        iterations+=1
        J = zeros(5,5)
        for k in 1:5
            step = zeros(5)
            step[k] = 1e-4
            J[:,k] = (last(at(v+step))-last(at(v-step)))/(2step[k])
        end
        d = -(pinv(J)*r)
        norm(d,Inf)>8 && (d .*= 8/norm(d,Inf))
        accepted = false
        for j in 0:20
            trial = v+0.5^j*d
            try
                gg,rr = at(trial)
                if norm(rr)<norm(r)
                    v,g,r = trial,gg,rr
                    accepted = true
                    break
                end
            catch err
                _recoverable(err) || rethrow()
            end
        end
        accepted || break
    end
    norm(r)<=2e-12 || throw(GeometryError("polygon solve stalled at residual $(norm(r)); try another initial tuple"))
    fine = geometry(g.e,2order)
    slit_error = maximum(abs.(fine.slit.-[h1,h2]))
    slit_error <= 1e-8 || throw(GeometryError("prescribed slit tips fail doubled-order check"))
    E = copy(g.e)
    for (j,level) in enumerate((m1,m2))
        trims[j]==0 && continue
        target = level + (choices[j]==:below ? -trims[j] : trims[j])
        E[_mouth(choices[j],j)] = _recover_end(g,_trim_wall(choices[j],j),target*T)
    end
    model = restore_model(e=g.e,endpoints=E,n=n,m1=m1,m2=m2,sides=choices,trims=trims,T=T,order=order)
    diagnostics = merge(model.diagnostics,(polygon_iterations=iterations,max_slit_error=slit_error))
    FilterModel(g,E,n,m1,m2,model.T,model.trims,choices,diagnostics)
end
