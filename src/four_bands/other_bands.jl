"""A marked elliptic Zolotarev solution whose alternation survives two passband gaps."""
struct ZolotarevFourBands
    endpoints::Vector{Float64}
    n::Int
    T::Float64
    q::Float64
    H::Float64
    K::Float64
    projective_scale::Float64
end

# Carlson's symmetric integral R_F. Duplication preserves its value and
# divides coordinate differences by four. At their mean the linear Taylor
# term vanishes by symmetry, hence a relative spread of 1e-8 suffices here.
function _zolo_rf(a,b,c)
    for _ in 1:64
        mean = (a+b+c)/3
        max(abs(a-mean),abs(b-mean),abs(c-mean)) <= 1e-8*abs(mean) &&
            return inv(sqrt(mean))
        sa,sb,sc = sqrt(a),sqrt(b),sqrt(c)
        λ = sa*(sb+sc)+sb*sc
        a,b,c = (a+λ)/4, (b+λ)/4, (c+λ)/4
    end
    throw(DomainError((a,b,c),"elliptic duplication failed to converge at this precision"))
end

function _zolo_normalize(f,x)
    A,D = first(f.endpoints),last(f.endpoints)
    A == -1 && D == 1 && f.projective_scale == 1 && return x
    left,right = f.projective_scale*(x-A), D-x
    (left-right)/(left+right)
end

function _zolo_physical(f,w)
    A,D = first(f.endpoints),last(f.endpoints)
    A == -1 && D == 1 && f.projective_scale == 1 && return w
    left,right = 1+w, f.projective_scale*(1-w)
    (D*left+A*right)/(left+right)
end

# On the wall w = q/dn(v K(q) | 1-q²), 0 ≤ v ≤ H.
function _zolo_tiles(f,w)
    t = clamp(abs(w),f.q,1.0)
    t==f.q && return 0.
    t==1 && return Float64(f.n)
    s = sqrt(((t-f.q)/t)*((t+f.q)/t)/((1-f.q)*(1+f.q)))
    c2 = ((1-t)/t)*((1+t)/t)*(f.q^2/((1-f.q)*(1+f.q)))
    s*_zolo_rf(c2,(f.q/t)^2,1.)/(f.K*f.T)
end

function _zolo_wall_x(f,j)
    j==0 && return f.endpoints[3]
    j==f.n && return f.endpoints[8]
    w = f.q*real(elliptic_value(1+im*(j*f.T),f.H))
    _zolo_physical(f,w)
end

function _zolo_model(e,n,q,H,scale)
    0<q<1 && isfinite(H) && H>0 && isfinite(scale) && scale>0 ||
        throw(ArgumentError("outer bands are unresolved at Float64 precision"))
    K = _zolo_rf(0.,(1-q)*(1+q),1.)
    f=ZolotarevFourBands(e,n,H/n,q,H,K,scale)
    _check_endpoints(e)
    # The first integer strictly after the lower gap end suffices. Inspect
    # its neighbours as well to allow for rounding when the end is an extremum.
    for i in (4,6)
        l,r = e[i],e[i+1]
        k = floor(Int,_zolo_tiles(f,_zolo_normalize(f,l)))
        for j in max(0,k-1):min(n,k+2)
            l<_zolo_wall_x(f,j)<r && throw(ArgumentError(
                "passband gap ($l,$r) removes the Zolotarev extremum at tile $j"))
        end
    end
    f
end

"""Synthesize the genus-one solution if both positive gaps preserve alternation."""
function zolotarev_fourbands(ends,n::Integer)
    n>0 || throw(ArgumentError("degree must be positive"))
    e=Float64.(collect(ends))
    _check_endpoints(e)
    A,B,C,D=e[[1,2,3,8]]
    s=sqrt((B-A)/(C-A))*sqrt((D-C)/(D-B))
    q=(((C-B)/(C-A))*((D-A)/(D-B)))/(1+s)^2
    0<q<1 && q^2>0 || throw(ArgumentError("unresolved two-band cross ratio"))
    H=_zolo_rf(0.,q^2,1.)/_zolo_rf(0.,(1-q)*(1+q),1.)
    _zolo_model(e,n,q,H,s*(D-B)/(B-A))
end

"""Construct canonical bands from τ and two removed open intervals in tile heights."""
function zolotarev_fourbands(tau::Complex,n::Integer; gaps=((.2,.4),(n-.4,n-.2)))
    iszero(real(tau)) && isfinite(imag(tau)) && imag(tau)>0 ||
        throw(ArgumentError("tau must be purely imaginary with positive imaginary part"))
    n>0 || throw(ArgumentError("degree must be positive"))
    length(gaps)==2 && all(g->length(g)==2,gaps) ||
        throw(ArgumentError("expected two pairs of tile heights"))
    j1,j2=gaps[1]
    j3,j4=gaps[2]
    0<j1<j2<j3<j4<n || throw(ArgumentError("gaps must lie in increasing order inside (0,n)"))
    for (l,r) in gaps
        floor(l)+1<r && throw(ArgumentError("open gap ($l,$r) removes an integer tile height"))
    end
    T=Float64(imag(tau))
    H=n*T
    q=elliptic_constants(H).k
    xs=[q*real(elliptic_value(1+im*(j*T),H)) for j in (j1,j2,j3,j4)]
    _zolo_model([-1.,-q,q,xs...,1.],n,q,H,1.)
end

function wall_values(f::ZolotarevFourBands,x::Real,b::Integer)
    1<=b<=4 || throw(ArgumentError("band index must be in 1:4"))
    f.endpoints[2b-1]<=x<=f.endpoints[2b] || throw(ArgumentError("point outside selected band"))
    _zolo_wall(f,x,b==1 ? -1 : 1)
end

function _zolo_wall(f,x,target)
    w = _zolo_normalize(f,x)
    e=f.endpoints
    tiles = x==e[1] || x==e[8] ? Float64(f.n) : x==e[2] || x==e[3] ? 0. : _zolo_tiles(f,w)
    c=elliptic_constants(f.T)
    phase=rem(tiles,2.,RoundNearest)
    ripple=_wall_ripple(phase,f.T,c,target)
    reciprocal_ripple=_wall_ripple(phase+1,f.T,c,target)
    a,b=_pair(target+im*phase*f.T,f.T,c)
    (; direct=real(a),reciprocal=real(b),ripple,reciprocal_ripple,target,height_tiles=tiles,method=:elliptic)
end

function evaluate(f::ZolotarevFourBands,z::Number;method=:elliptic)
    method===:elliptic || throw(ArgumentError("this family uses elliptic evaluation"))
    isfinite(z) || throw(ArgumentError("expected a finite evaluation point"))
    if isreal(z)
        x=real(z)
        e=f.endpoints
        e[1]<=x<=e[2] && return _zolo_wall(f,x,-1)
        e[3]<=x<=e[8] && return _zolo_wall(f,x,1)
        w=_zolo_normalize(f,x)
        top=abs(w)>1
        s=clamp(top ? inv(w) : w/f.q,-1.,1.)
        u=s*_zolo_rf((1-s)*(1+s),(1-f.q*s)*(1+f.q*s),1.)/f.K
        a,b=_pair(u,f.T,elliptic_constants(f.T))
        direct,reciprocal = top && isodd(f.n) ? (real(b),real(a)) : (real(a),real(b))
        return (;direct,reciprocal,method)
    end
    w=_zolo_normalize(f,complex(z))
    s=w/f.q
    u=s*_zolo_rf(1-s^2,1-w^2,1.)/f.K
    a,b=_pair(u,f.T,elliptic_constants(f.T))
    (;direct=a,reciprocal=b,method)
end

_plot_counts(f::ZolotarevFourBands) = ntuple(b -> b==1 ? f.n :
    ceil(Int,_zolo_tiles(f,_zolo_normalize(f,f.endpoints[2b])) -
             _zolo_tiles(f,_zolo_normalize(f,f.endpoints[2b-1]))),4)
_plot_description(::ZolotarevFourBands) = "Дробь Золотарёва · род 1"

"""Degree-three real Pick fraction with four active intervals and canonical genus five."""
struct PickFourBands
    endpoints::Vector{Float64}
    n::Int
    T::Float64
    poles::NTuple{2,Float64}
    residues::NTuple{2,Float64}
    offset::Float64
end

_pick_r(x,poles,residues,off) = x + off - residues[1]/(x-poles[1]) - residues[2]/(x-poles[2])

function _pick_preimages(v,poles,residues,t,mu)
    a,b=poles
    A,B=residues
    # The eigenvalue equation of this symmetric arrow matrix is r(x)=v.
    xs=eigvals(Symmetric([v-t sqrt(A) sqrt(B); sqrt(A) a 0.; sqrt(B) 0. b]))
    for i in eachindex(xs)
        x=xs[i]
        slope=1+A/(x-a)^2+B/(x-b)^2
        refined=x-(_pick_r(x,(a,b),(A,B),t)-v)/slope
        if abs(_pick_r(refined,(a,b),(A,B),t)-v)<abs(_pick_r(x,(a,b),(A,B),t)-v)
            xs[i]=refined
        end
    end
    maximum(x->abs(_pick_r(x,(a,b),(A,B),t)-v),xs)<=2e-8*mu ||
        throw(ArgumentError("band endpoints are unresolved relative to the requested ripple"))
    xs
end

"""
    pick_fourbands(tau; poles=(-1.,1.), residues=(1.,1.), offset=0.)

Take inverse images of the ripple levels under a strictly increasing real
rational cubic. The eight selected endpoints certify local optimality.
The remaining four odd preimages contribute to the source curve's genus.
"""
function pick_fourbands(tau::Complex; poles=(-1.,1.), residues=(1.,1.), offset=0.)
    iszero(real(tau)) && isfinite(imag(tau)) && imag(tau)>0 ||
        throw(ArgumentError("tau must be purely imaginary with positive imaginary part"))
    length(poles)==2 && length(residues)==2 ||
        throw(ArgumentError("expected two poles and two residues"))
    a,b=Float64.(poles)
    A,B=Float64.(residues)
    t=Float64(offset)
    all(isfinite,(a,b,A,B,t)) && a<b && A>0 && B>0 ||
        throw(ArgumentError("expected two increasing finite poles, positive residues and a finite offset"))
    T=Float64(imag(tau))
    mu=elliptic_constants(T).mu
    0<mu<1 || throw(ArgumentError("four distinct ripple levels must be resolved in Float64"))
    preimages(v)=_pick_preimages(v,(a,b),(A,B),t,mu)
    lo,hi=preimages(1-mu),preimages(1+mu)
    e=[preimages(-1-mu)[1],preimages(-1+mu)[1],lo[1],hi[1],lo[2],hi[2],lo[3],hi[3]]
    _check_endpoints(e)
    PickFourBands(e,3,T,(a,b),(A,B),t)
end

"""Sorted twelve simple branch points of the canonical double cover."""
function pick_branchpoints(f::PickFourBands)
    mu=elliptic_constants(f.T).mu
    sort!(vcat((_pick_preimages(v,f.poles,f.residues,f.offset,mu)
                for v in (-1-mu,-1+mu,1-mu,1+mu))...))
end

function evaluate(f::PickFourBands,z::Number;method=:rational)
    method===:rational || throw(ArgumentError("this family uses explicit rational evaluation"))
    isfinite(z) || throw(ArgumentError("expected a finite evaluation point"))
    z in f.poles && throw(DomainError(z,"direct pole"))
    direct=_pick_r(z,f.poles,f.residues,f.offset)
    iszero(direct) && throw(DomainError(z,"reciprocal pole"))
    mu=elliptic_constants(f.T).mu
    reciprocal=((1-mu)*(1+mu))/direct
    (;direct,reciprocal,method)
end

_plot_counts(::PickFourBands) = (1,1,1,1)
_plot_description(::PickFourBands) = "Рациональная дробь степени 3 · род 5"

"""A four-band cubic with a simple tangency at the lower stop level, of genus four."""
struct FoldedFourBands
    endpoints::Vector{Float64}
    n::Int
    T::Float64
    cut::Float64
    branchpoints::Vector{Float64}
    scale::Float64
    turning::Float64
    odd_root::Float64
end

# There is one root in each bracket. For x > 2*turning the rational
# term (x-odd_root)*(x-turning)²/(x²-1) exceeds x/4.
function _folded_roots(delta,scale,turning,odd_root)
    level=delta/scale
    upper=max(2turning,4level)+1
    isfinite(upper) || throw(ArgumentError("cubic branch bounds overflow Float64"))
    polynomial(x)=(x-odd_root)*(x-turning)^2-level*(x-1)*(x+1)
    map(((odd_root,-1.),(1.,turning),(turning,upper))) do (a,b)
        signleft=signbit(polynomial(a))
        for _ in 1:64
            c=(a+b)/2
            fc=polynomial(c)
            iszero(fc) && return c
            (c==a || c==b) && return abs(polynomial(a))<abs(polynomial(b)) ? a : b
            if signbit(fc)==signleft
                a=c
            else
                b=c
            end
        end
        (a+b)/2
    end
end

"""Construct the genus-four cubic from τ using one stop tangency and three full passbands."""
function folded_fourbands(tau::Complex; scale=1.,turning=2.,odd_root=-3.)
    iszero(real(tau)) && isfinite(imag(tau)) && imag(tau)>0 ||
        throw(ArgumentError("tau must be purely imaginary with positive imaginary part"))
    T=Float64(imag(tau))
    mu=elliptic_constants(T).mu
    0<mu<1 || throw(ArgumentError("four distinct ripple levels must be resolved in Float64"))
    all(x->x isa Real && !(x isa Bool) && isfinite(x),(scale,turning,odd_root)) &&
        scale>0 && turning>1 && odd_root < -1 ||
        throw(ArgumentError("expected scale > 0, turning > 1, odd_root < −1"))
    a,t,d=Float64(scale),Float64(turning),Float64(odd_root)
    roots(delta)=_folded_roots(delta,a,t,d)
    stop,lo,hi=roots(2mu),roots(2.),roots(2+2mu)
    cut=roots(1+mu)[2] # R=0 between the last passband and the stopband
    transform(x)=-inv(x-cut)
    e=transform.([stop[2],stop[3],lo[3],hi[3],lo[1],hi[1],hi[2],lo[2]])
    _check_endpoints(e)
    branchpoints=sort(transform.([d,stop...,lo...,hi...]))
    f=FoldedFourBands(e,3,T,cut,branchpoints,a,t,d)
    for b in 1:4, j in 2b-1:2b
        x=e[j]
        v=evaluate(f,x).direct
        abs(abs(v-(b==1 ? -1 : 1))-mu)<=2e-8*mu ||
            throw(ArgumentError("band endpoints are unresolved relative to the requested ripple"))
    end
    f
end

function evaluate(f::FoldedFourBands,z::Number;method=:rational)
    method===:rational || throw(ArgumentError("this family uses explicit rational evaluation"))
    isfinite(z) || throw(ArgumentError("expected a finite evaluation point"))
    iszero(z) && throw(DomainError(z,"direct pole"))
    x=f.cut-inv(z)
    abs(x)==1 && isreal(x) && throw(DomainError(z,"direct pole"))
    mu=elliptic_constants(f.T).mu
    direct=f.scale*(x-f.odd_root)*(x-f.turning)^2/((x-1)*(x+1))-(1+mu)
    iszero(direct) && throw(DomainError(z,"reciprocal pole"))
    reciprocal=((1-mu)*(1+mu))/direct
    (;direct,reciprocal,method)
end

function wall_values(f::Union{FoldedFourBands,PickFourBands},x::Real,b::Integer)
    1<=b<=4 || throw(ArgumentError("band index must be in 1:4"))
    f.endpoints[2b-1]<=x<=f.endpoints[2b] || throw(ArgumentError("point outside selected band"))
    v=evaluate(f,x)
    mu=elliptic_constants(f.T).mu
    target=b==1 ? -1 : 1
    ripple=(v.direct-target)/mu
    reciprocal_ripple=(-mu-target*ripple)/v.direct
    (;v...,target,ripple,reciprocal_ripple)
end

_plot_counts(::FoldedFourBands) = (2,1,1,1)
_plot_description(::FoldedFourBands) = "Рациональная дробь степени 3 · род 4"
