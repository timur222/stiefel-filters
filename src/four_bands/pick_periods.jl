# Periods of y²=∏(x−eᵣ) for ordered marked real branch points.
# The even intervals are the a-ovals. The b-periods are cumulative signed
# sums along odd intervals in the same marking as the genus-three module.
function _branch_omega(e,order)
    nb=length(e)
    iseven(nb) && 4<=nb<=12 && all(isfinite,e) && all(>(0),diff(e)) ||
        throw(ArgumentError("expected 2g+2 distinct ordered real branch points, 1 ≤ g ≤ 5"))
    g=(nb-2)÷2
    nodes,weights=_gauss(order)
    a,b=first(e),last(e)
    v=(e .- (a/2+b/2)) ./ (b/2-a/2)
    all(>(0),diff(v)) || throw(GeometryError("unresolved branch spacing in the period coordinate"))
    M=zeros(nb-1,g)
    for s in 1:nb-1, k in eachindex(nodes)
        # x=a+(b−a)sin²(t), t∈[0,π/2]. Compute small distances from
        # the closer endpoint directly to avoid cancellation near each root.
        angle=(π/4)*(1+nodes[k])
        left=(v[s+1]-v[s])*sin(angle)^2
        right=(v[s+1]-v[s])*cos(angle)^2
        y=1.
        for l in 1:s-1
            y *= sqrt(v[s]-v[l]+left)
        end
        for l in s+2:nb
            y *= sqrt(v[l]-v[s+1]+right)
        end
        x=v[s]+left
        w=(π/2)*weights[k]/y
        M[s,1]+=w
        g>=2 && (M[s,2]+=w*x)
        previous,current=1.,x
        for j in 3:g
            previous,current=current,2x*current-previous
            M[s,j]+=w*current
        end
    end
    signs=[isodd(i) ? -1. : 1. for i in 1:g]
    A=M[2:2:2g,:] .* signs
    B=cumsum(M[1:2:2g-1,:] .* signs;dims=1)
    Y=B/A
    discrepancy=maximum(abs,Y-transpose(Y))
    all(isfinite,Y) && discrepancy <= 2e-8*max(1.,maximum(abs,Y)) ||
        throw(GeometryError("period integration lost symmetry; increase order or check branch spacing"))
    Ω=(Y+transpose(Y))/2
    isposdef(Symmetric(Ω)) || throw(GeometryError("period integration lost positivity"))
    return Ω
end

"""
    pick_period_matrix(f::PickFourBands; order=128)

Real imaginary part Ω of the marked period matrix Π=iΩ of the cubic Pick
pullback. All twelve preimages of the four ripple levels enter the integral.
The marking uses the ascending branch order and five finite a-ovals on
[e₂,e₃],…,[e₁₀,e₁₁].
"""
function pick_period_matrix(f::PickFourBands; order::Int=128)
    _branch_omega(pick_branchpoints(f),order)
end

function _marked_target(P,g::Int)
    size(P)==(g,g) && all(v->v isa Number && isfinite(v),P) ||
        throw(ArgumentError("expected a finite marked $(g)×$(g) Ω or Π=iΩ"))
    Ω = if all(isreal,P)
        Float64.(real.(P))
    elseif all(x->iszero(real(x)),P)
        Float64.(imag.(P))
    else
        throw(ArgumentError("expected real Ω or purely imaginary Π=iΩ"))
    end
    all(isfinite,Ω) && maximum(abs,Ω-transpose(Ω))<=5e-13*max(1.,maximum(abs,Ω)) ||
        throw(ArgumentError("expected finite symmetric Ω"))
    Ω=(Ω+transpose(Ω))/2
    isposdef(Symmetric(Ω)) || throw(ArgumentError("expected positive definite Ω"))
    return Ω
end

function _pick_seed(p)
    isnothing(p) && return zeros(4)
    if p isa PickFourBands
        a,b=p.poles
        p=((b-a)/2,p.residues...,p.offset+(a+b)/2)
    end
    length(p)==4 && all(v->v isa Real && isfinite(v),p) && all(>(0),(p[1],p[2],p[3])) ||
        throw(ArgumentError("initial must give positive pole halfspacing, two positive residues, and offset"))
    return [log(p[1]),log(p[2]),log(p[3]),p[4]]
end

# Damped Gauss-Newton in real coordinates with one residual per marked period.
function _fit_marked_periods(at,v)
    f,r=at(v)
    for iter in 1:70
        norm(r,Inf)<1e-11 && break
        J=zeros(length(r),length(v))
        for i in eachindex(v)
            h=1e-4*max(1.,abs(v[i]))
            dv=zeros(length(v))
            dv[i]=h
            J[:,i]=(last(at(v+dv))-last(at(v-dv)))/(2h)
        end
        step=-pinv(J)*r
        all(isfinite,step) || break
        norm(step,Inf)>2 && (step .*= 2/norm(step,Inf))
        advanced=false
        for k in 0:16
            next=v+exp2(-k)*step
            try
                ff,rr=at(next)
                if norm(rr)<norm(r)
                    v,f,r=next,ff,rr
                    advanced=true
                    break
                end
            catch ex
                ex isa GeometryError || ex isa ArgumentError || ex isa DomainError || rethrow()
            end
        end
        advanced || break
    end
    return f,r
end

"""
    pick_from_periods(tau, omega; order=128, initial=nothing)

Local inverse in the marked degree-three Pick family. Normalize the poles to
(−s,s) and the slope at infinity to 1. Fit s, both positive residues and the
offset to all fifteen independent entries. Re-integrate at twice the order;
return a filter only if both matrix residuals are ≤1e−8 when scaled by
max(|Ωᵢⱼ|,0.1). An optional seed is a model or (s,residue₁,residue₂,offset).
The returned coordinate may differ projectively from the original bands.
"""
function pick_from_periods(tau,P::AbstractMatrix; order::Int=128,initial=nothing)
    T=_period_tau(tau)
    target=_marked_target(P,5)
    # Pull back dz / √∏(z−vⱼ) from the elliptic target. The five a-cycles
    # project with alternating signs. On the b-cycles the cumulative two
    # adjacent cuts cancel, leaving one elliptic period at each odd index.
    s=[1.,-1.,1.,-1.,1.]
    expected=(T/2).*[1.,0.,1.,0.,1.]
    norm(target*s-expected,Inf)<=1e-8*max(1.,T,norm(target,Inf)) ||
        throw(ArgumentError("periods incompatible with the cubic elliptic covering and τ"))
    index=[CartesianIndex(i,j) for i in 1:5 for j in i:5]
    scale=max.(.1,abs.(target[index]))
    order>=2 || throw(ArgumentError("quadrature order must be at least two"))
    function at(v)
        s,c,d=exp.(v[1:3])
        f=pick_fourbands(im*T;poles=(-s,s),residues=(c,d),offset=v[4])
        return f,(pick_period_matrix(f;order)[index]-target[index])./scale
    end
    f,r=_fit_marked_periods(at,_pick_seed(initial))
    fine=pick_period_matrix(f;order=2order)
    residual=max(norm(r,Inf),maximum(abs,(fine[index]-target[index])./scale))
    residual<=1e-8 || throw(GeometryError("genus-five Pick reconstruction failed (scaled full-matrix residual $residual); check the marking, family, τ, seed and integration order"))
    return f
end
