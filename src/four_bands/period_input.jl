"""Period data in the marked basis of `ThetaCurve`; τ is purely imaginary and Π=iΩ."""
function period_parameters(f::FilterModel)
    omega = ThetaCurve(f.geom).omega
    (; tau=im*f.T, omega, n=f.n, m1=f.m1, m2=f.m2,
       omega22=omega[2,2], omega23=omega[2,3], omega33=omega[3,3],
       trims=f.trims, sides=f.sides)
end

function _period_tau(tau)
    tau isa Number && !(tau isa Bool) && isfinite(tau) ||
        throw(ArgumentError("expected a finite τ on the positive imaginary axis"))
    T = tau isa Real ? tau : (iszero(real(tau)) ? imag(tau) : NaN)
    isfinite(Float64(T)) && T>0 || throw(ArgumentError("expected τ=iT with T>0 (or the real number T)"))
    Float64(T)
end

function _marked_omega(matrix)
    size(matrix)==(3,3) || throw(ArgumentError("the two-slit family needs a marked 3×3 matrix Ω or Π=iΩ"))
    all(v->v isa Number && isfinite(v),matrix) || throw(ArgumentError("nonfinite period entry"))
    omega = if all(isreal,matrix)
        Float64.(real.(matrix))
    elseif all(v->iszero(real(v)),matrix)
        Float64.(imag.(matrix))
    else
        throw(ArgumentError("in this marking Π must be purely imaginary (or supply real Ω)"))
    end
    all(isfinite,omega) || throw(ArgumentError("period matrix outside Float64 range"))
    maximum(abs,omega-transpose(omega)) <= 5e-13*max(1.,maximum(abs,omega)) ||
        throw(ArgumentError("Ω must be symmetric"))
    omega = (omega+transpose(omega))/2
    isposdef(Symmetric(omega)) || throw(ArgumentError("Ω must be positive definite"))
    omega
end

function _period_mode(T,omega)
    q = (2/T).*omega[1,:]
    all(x->isfinite(x) && 0<x<2.0^53,q) || throw(ArgumentError("first row gives inadmissible integer heights"))
    k = round.(Int,q)
    maximum(abs,q.-k)<=1e-7 || throw(ArgumentError("first row is incompatible with τ: require 2Ω₁ⱼ / Im(τ) integral"))
    n = k[1]
    (n,n-k[2],n-k[3])
end

# Recover ordered branches in the projective gauge (-1,0,1). In the
# five-equation version Ω₂₃ is an OUTPUT. In the six-equation version its
# residual supplies the independent condition on the hyperelliptic locus.
function _fit_marked_omega(target,known,order,initial)
    scales = max.(abs.(target[known]),0.1)
    function at(v)
        g = geometry(_polygon_points(v),order)
        o = ThetaCurve(g).omega
        g,(o[known]-target[known])./scales
    end
    v = _polygon_seed(initial)
    g,r = at(v)
    iterations = 0
    while norm(r,Inf)>5e-14 && iterations<80
        iterations += 1
        J = zeros(length(r),5)
        for k in eachindex(v)
            d = zeros(5)
            d[k] = 1e-4
            J[:,k]=(last(at(v+d))-last(at(v-d)))/(2d[k])
        end
        step = -(pinv(J)*r)
        all(isfinite,step) || break
        size_step = norm(step,Inf)
        size_step>5 && (step .*= 5/size_step)
        accepted = false
        for j in 0:18
            trial = v + 0.5^j*step
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
    # Independent period integration at twice the quadrature order.
    fine = ThetaCurve(geometry(g.e,2order)).omega
    residual = maximum(abs,(fine[known]-target[known])./scales)
    max(norm(r,Inf),residual)<=1e-9 ||
        throw(GeometryError("period reconstruction failed (scaled residual $(max(norm(r,Inf),residual))); check the marked hyperelliptic locus, initialization and quadrature order"))
    g,iterations,residual, fine
end

function _from_period_data(T,mode,omega,has_cross,trims,sides,order,initial)
    n,m1,m2 = mode
    all(v->v isa Integer && !(v isa Bool),mode) && 0<m1<m2<n<2^53 ||
        throw(ArgumentError("require integer 0<m1<m2<n<2^53"))
    length(trims)==2 && all(v->v isa Real && isfinite(v) && 0<=v<=1,trims) ||
        throw(ArgumentError("expected two trim lengths in [0,1]"))
    length(sides)==2 && all(s->s in (:above,:below,"above","below"),sides) ||
        throw(ArgumentError("expected two trim sides, above or below"))
    order>=2 || throw(ArgumentError("quadrature order must be at least two"))
    known = CartesianIndex.([(1,1),(1,2),(1,3),(2,2),(3,3)])
    has_cross && push!(known,CartesianIndex(2,3))
    g,iterations,residual,fine = _fit_marked_omega(omega,known,order,initial)
    choices = Symbol.(Tuple(sides))
    E = copy(g.e)
    for (j,level) in enumerate((m1,m2))
        iszero(trims[j]) && continue
        height = level + (choices[j]==:below ? -trims[j] : trims[j])
        E[_mouth(choices[j],j)] = _recover_end(g,_trim_wall(choices[j],j),height*T)
    end
    # This checks all four wall heights in units of the requested τ and
    # validates slit positions and retained bands, again at doubled order.
    f = restore_model(e=g.e,endpoints=E,n=Int(n),m1=Int(m1),m2=Int(m2),
                      sides=choices,trims=trims,T=T,order=order)
    checks = merge(f.diagnostics,(period_iterations=iterations,period_scaled_residual=residual,
                                  omega23=fine[2,3]))
    FilterModel(f.geom,f.endpoints,f.n,f.m1,f.m2,f.T,f.trims,f.sides,checks)
end

"""
    filter_from_periods(tau, omega; trims=(0.,0.), sides=(:above,:above), order=128, initial=nothing)

Invert real Ω (or purely imaginary Π=iΩ) in the marked genus-three basis.
The first row determines integer n,m₁,m₂ via 2Ω₁⋅=Im(τ)(n,n−m₁,n−m₂).
Recover eight ordered real branch points with e₁=−1, e₂=0, e₈=1, and check
all six measured entries after doubling the integration order. Positive real
`tau` denotes Im(τ). `initial` is an optional ordered set of eight seed points.
A failed local inversion is reported as an error.
"""
function _three_from_periods(tau,matrix::AbstractMatrix;trims=(0.,0.),sides=(:above,:above),
                             order::Int=128,initial=nothing)
    T = _period_tau(tau)
    omega = _marked_omega(matrix)
    mode = _period_mode(T,omega)
    _from_period_data(T,mode,omega,true,trims,sides,order,initial)
end

"""
    filter_from_periods(tau,n,m1,m2; omega22,omega33,omega23=nothing,kwargs...)

Build the quantized first row from the degree and slit levels. When Ω₂₃ is
omitted, solve five period equations in the five real branch moduli and
read Ω₂₃ from the resulting curve. The result is a local numerical solution.
If Ω₂₃ is given, all six equations must be consistent to scaled tolerance
1e−9 (scale max(0.1,|Ωᵢⱼ|)). Wall quantization is checked to 1e−7 tiles.
"""
function filter_from_periods(tau,n::Integer,m1::Integer,m2::Integer;omega22,omega33,omega23=nothing,
                             trims=(0.,0.),sides=(:above,:above),order::Int=128,initial=nothing)
    T = _period_tau(tau)
    all(v->!(v isa Bool),(n,m1,m2)) && 0<m1<m2<n<2^53 ||
        throw(ArgumentError("require integer 0<m1<m2<n<2^53"))
    all(v->v isa Real && isfinite(v) && v>0,(omega22,omega33)) ||
        throw(ArgumentError("Ω₂₂, Ω₃₃ must be finite and positive"))
    cross = isnothing(omega23) ? 0. : omega23
    cross isa Real && isfinite(cross) || throw(ArgumentError("Ω₂₃ must be real and finite"))
    a,b,c = (T/2).*(n,n-m1,n-m2)
    omega = Float64[a b c; b omega22 cross; c cross omega33]
    all(isfinite,omega) || throw(ArgumentError("period matrix outside Float64 range"))
    # Schur complement permits some positive-definite completion exactly
    # when both of its diagonal entries are positive.
    b^2<a*omega22 && c^2<a*omega33 || throw(ArgumentError("principal minors of Ω must be positive"))
    !isnothing(omega23) && _marked_omega(omega)
    _from_period_data(T,(n,m1,m2),omega,!isnothing(omega23),trims,sides,order,initial)
end
