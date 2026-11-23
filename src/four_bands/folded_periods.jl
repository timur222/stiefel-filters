"""
    folded_period_matrix(f; order=128)

Imaginary part Ω of the marked genus-four period matrix Π=iΩ. The four
real a-cycles lie above the even intervals in the ascending physical
branch order. The tangent double preimage is absent from this list.
"""
folded_period_matrix(f::FoldedFourBands;order::Int=128) = _branch_omega(f.branchpoints,order)

"""
    folded_from_periods(tau, omega; order=128, initial=nothing)

Fit the parameters of the marked genus-four cubic to all ten entries, reducing
the step when needed. Recompute Ω at twice the quadrature order and accept
only if the scaled maximum residual is at most 1e-8 for both computations.
Optionally seed with a model or with (scale, turning, odd_root).
"""
function folded_from_periods(tau,P::AbstractMatrix;order::Int=128,initial=nothing)
    T=_period_tau(tau)
    target=_marked_target(P,4)
    order>=2 || throw(ArgumentError("quadrature order must be at least two"))
    index=[CartesianIndex(i,j) for i in 1:4 for j in i:4]
    scale=max.(.1,abs.(target[index]))
    function at(v)
        a,t,d=exp.(v)
        f=folded_fourbands(im*T;scale=a,turning=1+t,odd_root=-1-d)
        f,(folded_period_matrix(f;order)[index]-target[index])./scale
    end
    seed = if isnothing(initial)
        (1.,2.,-3.)
    elseif initial isa FoldedFourBands
        (initial.scale,initial.turning,initial.odd_root)
    else
        initial
    end
    length(seed)==3 && all(x->x isa Real && isfinite(x),seed) &&
        seed[1]>0 && seed[2]>1 && seed[3] < -1 ||
        throw(ArgumentError("expected initial (scale > 0, turning > 1, odd_root < −1)"))
    v=log.([seed[1],seed[2]-1,-1-seed[3]])
    f,r=_fit_marked_periods(at,v)
    fine=folded_period_matrix(f;order=2order)
    residual=max(norm(r,Inf),maximum(abs,(fine[index]-target[index])./scale))
    residual<=1e-8 || throw(GeometryError("genus-four cubic reconstruction failed (scaled full-matrix residual $residual); check marking, family, τ, seed and quadrature order"))
    f
end
