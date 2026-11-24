# Two rectangular charts of unit width. The tile heights in these
# coordinates are T_i/2. The phase at the shared edge is m+offset mod 2.
struct WeldBlock
    widthfactor::Float64
    delta::Float64
    k::Float64
    K::Float64
    T::Float64
    c::Float64
    edge::Matrix{Float64}
    seed::NamedTuple
end

# u is deviation in units of the symmetric two-band Zolotarev deviation.
# An increasing real Möbius transformation maps its four extremal values
# to -1±ε and 1±bε. On each band it acts on u via (u+h)/(1+h*u).
function _weighted_phase(c,T,seed,phase,target)
    u=_wall_ripple(rem(phase,2.,RoundNearest),T,seed,target)
    h=c*seed.mu/(1+target*c)
    (u+h)/(1+h*u)
end

function _weld_block(e,b,odd)
    d=2sqrt(b)*e/sqrt(4-((b-1)*e)^2)
    0<d<1 || throw(ArgumentError("the two ripple intervals must be strictly disjoint"))
    k=(1-d)/(1+d)
    0<k<1 || throw(ArgumentError("the ripple levels coalesce at Float64 precision"))
    K=_zolo_rf(0.,4d/(1+d)^2,1.)
    T=_zolo_rf(0.,k^2,1.)/K
    c=-(b-1)/(b+1+hypot(2sqrt(b),(b-1)*d))
    A=c+b*e*(1+c)/d
    F=[A 1+c+b*e*c*d-A; c 1.]
    odd && (F=F*[0. 1-d^2; 1. 0.])
    WeldBlock(b,d,k,K,T,c,F,elliptic_constants(T))
end

# eta parametrizes the common edge in the narrow block, using the entire
# edge of the broad block. Both coordinates increase from the stop wall.
function _weld_edge(b,a,A,s)
    if s<1e-5
        den=A[2,2]-A[2,1]*(1-b.delta)
        slope=(b.K/a.K)*sqrt(det(A)*(1-b.delta)*(1-b.k^2)/
                    (den^2*(1-a.delta)*(1-a.k^2)))
        return (s*slope,slope)
    end
    sn=s==1 ? 1. : real(elliptic_value(2s-1,b.T))
    z=(1-b.delta)*sn
    den=A[2,1]*z+A[2,2]
    r=(A[1,1]*z+A[1,2])/(den*(1-a.delta))
    abs(r)<=1+1e-11 || throw(DomainError(r,"incompatible seam coordinates"))
    r=clamp(r,-1.,1.)
    r2=(1-r)*(1+r); dn2=(1-a.k*r)*(1+a.k*r)
    eta=abs(r)==1 ? (1+r)/2 : (1+r*_zolo_rf(r2,dn2,1.)/a.K)/2
    cn=sqrt(max(0.,(1-sn)*(1+sn)))
    der=s==1 ? 0. : (1-b.delta)*b.K*cn*sqrt((1-b.k*sn)*(1+b.k*sn))*det(A)/
        (den^2*(1-a.delta)*a.K*sqrt(r2*dn2))
    eta,der
end

"""
    welded_threebands(n,m,epsilon; third_weight=2, offset=0, order=24)

Direct welded quadrilateral in the mode m, with n−m tiles in block 2.
The field modulus is the height divided by the width of the uniformized
rectangle. Fourier order controls a one-dimensional interface calculation.
"""
function welded_threebands(n::Integer,m::Integer,e::Real;third_weight=2.,offset=0,order::Int=24)
    1<=m<n<2^52 && offset in (0,1) && order>=4 || throw(ArgumentError("expected 0<m<n<2^52, offset 0 or 1, order ≥4"))
    rho=1/Float64(third_weight); e=Float64(e)
    isfinite(rho) && rho>0 && rho!=1 && 0<e<2/(1+max(1.,rho)) ||
        throw(ArgumentError("expected unequal positive widths with separate positive and negative bands"))
    blocks=(_weld_block(e,1.,isodd(m+offset)),_weld_block(e,rho,isodd(m+offset)))
    h=(m*blocks[1].T/2,(n-m)*blocks[2].T/2)
    bi=rho<1 ? 1 : 2; ai=3-bi
    b,a=blocks[bi],blocks[ai]; hb,ha=h[bi],h[ai]
    trans=a.edge\b.edge
    # The narrow block's infinite-strip logarithmic kernel minus the
    # identity-welding kernel is smooth even along the diagonal.
    count=max(48,3order,ceil(Int,24/(π*ha)))
    s=[(j-.5)/count for j in 1:count]
    edge=[_weld_edge(b,a,trans,v) for v in s]
    eta=first.(edge); der=last.(edge)
    z=cospi.(s); v=cospi.(eta)
    kernel=[-log(i==j ? sinpi(eta[i])*der[i]/sinpi(s[i]) :
                (v[i]-v[j])/(z[i]-z[j]))/π for i in 1:count,j in 1:count]
    D=[cospi(j*t) for j in 0:order,t in s]
    S=D*kernel*transpose(D)/count^2
    for j in 1:order
        S[j+1,j+1]+=(1+tanh(π*j*hb))/(2π*j)
    end
    ks=1:max(1,ceil(Int,22/(π*ha)))
    C=D*[cospi(k*v) for v in eta,k in ks]/count
    S .+= (2/π).*(C*Diagonal([(tanh(π*k*ha)-1)/k for k in ks])*transpose(C))
    coeff=-(cholesky(Symmetric(S[2:end,2:end]))\S[2:end,1])
    excess=S[1,1]+dot(S[1,2:end],coeff)
    modulus=(hb+ha)+excess
    isfinite(modulus) && modulus>0 || throw(ArgumentError("unresolved welded modulus"))
    (;n,m,epsilon=e,third_weight,offset,order,blocks,h,bi,ai,trans,coeff,modulus,excess,
        edge_end=first(_weld_edge(b,a,trans,1.)))
end

function _weld_sinh_ratio(k,y,h)
    exp(-π*k*(h-y))*(-expm1(-2π*k*y))/(1+exp(-2π*k*h))
end

# Cache the narrow-block seam quadrature for repeated wall evaluations.
# Cubic changes of variable resolve logarithms at either end of the seam.
function _weld_grid(f,side)
    b,a=f.blocks[f.bi],f.blocks[f.ai]
    nodes,weights=_gauss(max(96,4f.order))
    t=(nodes .+ 1)./2
    u=t.^3
    s=side==:left ? u : side==:tip ? 1 .- u : t
    jac=side==:right ? ones(length(t)) : 3 .* t.^2
    eta=[first(_weld_edge(b,a,f.trans,v)) for v in s]
    pd=[1+sum(f.coeff[j]*cospi(j*v) for j in 1:f.order) for v in s]
    mass=(weights./2).*jac.*pd
    ks=1:max(1,ceil(Int,34/(π*f.h[f.ai])))
    moments=[sum(mass.*cospi.(k.*eta)) for k in ks]
    # edge_end-eta(1-u) = curvature*u² + O(u⁴). Use this distance
    # directly whenever two rounded values of eta would nearly coincide.
    distance=f.edge_end .- eta
    if side==:tip
        A=f.trans
        den=A[2,1]*(1-b.delta)+A[2,2]
        limit=((A[1,1]*(1-b.delta)+A[1,2])/den)/(1-a.delta)
        curvature=(1-b.delta)*b.K^2*(1-b.k^2)*det(A)/
            (den^2*(1-a.delta)*a.K*sqrt((1-limit^2)*(1-(a.k*limit)^2)))
        for j in eachindex(u)
            u[j]<1e-4 && (distance[j]=curvature*u[j]^2)
        end
    end
    (;s,eta,mass,moments,distance)
end

# Harmonic conjugate measured from the remote horizontal side.
function _weld_remote(f,i,s0,y)
    h=f.h[i]; y==0 && return 0.
    exp(-π*(h-y)) == 0 && return y
    if i==f.bi
        return y+sum(f.coeff[j]*cospi(j*s0)*_weld_sinh_ratio(j,y,h)/(π*j) for j in 1:f.order)
    end
    side=s0==0 ? :left : y==h ? :tip : :right
    g=hasproperty(f,:grids) ? getproperty(f.grids,side) : _weld_grid(f,side)
    r=exp(-π*(h-y))
    corr=sum((_weld_sinh_ratio(k,y,h)-r^k)*cospi(k*s0)*g.moments[k]/k for k in eachindex(g.moments))
    dist(v)=hypot(1-r,2sqrt(r)*sinpi(v/2))
    for j in eachindex(g.mass)
        gap=side==:tip ? s0-f.edge_end+g.distance[j] : s0-g.eta[j]
        logs=-(log(dist(gap))+log(dist(s0+g.eta[j])))/2
        s0==0 && (logs+=log(dist(g.s[j])))
        corr+=g.mass[j]*logs
    end
    y+2corr/π+(s0==0 ? sum(f.coeff[j]*r^j/(π*j) for j in 1:f.order) : 0.)
end

"""Height in the uniformized rectangle at tile index t (counted from the bottom)."""
function welded_wall_height(f,i,wall,t)
    wall in (-1,1) && i in (1,2) || throw(ArgumentError("expected wall ±1 and block 1 or 2"))
    (i==1 ? 0<=t<=f.m : f.m<=t<=f.n) || throw(DomainError(t,"tile outside the requested block"))
    s0=(wall+1)/2
    if i==f.ai && wall== -1 && t==f.m
        return welded_wall_height(f,f.bi,wall,t)
    end
    y=(i==1 ? t : f.n-t)*f.blocks[i].T/2
    z=_weld_remote(f,i,s0,y)
    i==1 ? z : f.modulus-z
end

# q is the elliptic modulus, not the nome. The original Carlson inversion
# uses q². Switch coordinates before that square becomes subnormal.
function _weld_outer(M)
    if log(4.0) - π*M >= log(sqrt(floatmin(Float64)))
        outer = elliptic_constants(2M)
        return (; M, q=outer.k, K=outer.K)
    end
    q = 4 * exp(-big(π) * M)
    (; M, q, K=π/2)
end

# X(v)X(M-v)=q on the positive wall. For v≤M/2 the theta series gives
# X(v)=q*cosh(πv)*(1+O(exp(-πM))). At the precision switch the omitted
# relative correction is smaller than 1e-153, throughout this half-wall.
function _weld_abscissa(outer, v)
    M, q = outer.M, outer.q
    v == 0 && return q
    v == M && return one(q)
    if q isa BigFloat
        return v <= M/2 ? q * cosh(big(π)*v) : inv(cosh(big(π)*(M-v)))
    end
    q * real(elliptic_value(1+2im*v, 2M))
end

function _weld_height(outer, x)
    M, q = outer.M, outer.q
    a = abs(x)
    a == q && return 0.0
    a == 1 && return M
    if q isa BigFloat
        a = BigFloat(a)
        height = a <= sqrt(q) ? acosh(a/q)/big(π) : M - acosh(inv(a))/big(π)
        return clamp(Float64(height), 0.0, M)
    end
    tiles = _zolo_tiles((; n=2M, q, K=outer.K, T=1.0), a)
    clamp(Float64(tiles)/2, 0.0, M)
end

"""Parametric (x,r) on a vertical wall, in the symmetric outer-endpoint gauge."""
function welded_wall(f,i,wall,t,outer=_weld_outer(f.modulus))
    v=welded_wall_height(f,i,wall,t)
    x=wall*_weld_abscissa(outer,v)
    b=f.blocks[i]
    z=_weighted_phase(b.c,b.T,b.seed,t+f.offset,wall)
    (;x,value=wall+f.epsilon*(wall==1 ? b.widthfactor : 1.)*z)
end

function _weld_bisect(fun,goal,lo,hi,increasing)
    for _ in 1:60
        mid=(lo+hi)/2
        (mid==lo || mid==hi) && return mid
        level=fun(mid)
        level==goal && return mid
        if (level<goal)==increasing
            lo=mid
        else
            hi=mid
        end
    end
    (lo+hi)/2
end

# In a straight block the conjugate height is the tile height plus a
# bounded seam correction. Bracket that correction around the straight
# estimate, so the precision of an interior sample does not depend on how
# many tiles lie between it and the outer end.
function _weld_tile_at_height(f, block, wall, v)
    tile_height = f.blocks[block].T / 2
    if block == 1
        lo, hi = 0.0, Float64(f.m)
        guess = v / tile_height
    else
        lo, hi = Float64(f.m), Float64(f.n)
        guess = f.n - (f.modulus - v) / tile_height
    end
    t = clamp(guess, lo, hi)
    height(t) = welded_wall_height(f, block, wall, t)
    level = height(t)
    level == v && return t
    step = 1.0
    if level < v
        right = min(hi, t + step)
        while right < hi && height(right) < v
            step *= 2
            right = min(hi, t + step)
        end
        return _weld_bisect(height, v, t, right, true)
    end
    left = max(lo, t - step)
    while left > lo && height(left) > v
        step *= 2
        left = max(lo, t - step)
    end
    _weld_bisect(height, v, left, t, true)
end

"""Evaluate the welded rational function on a prescribed working interval."""
function welded_value(solution,x::Real)
    f=solution.filter; E=solution.ends; target=solution.target
    band=findfirst(j->E[2j-1]<=x<=E[2j],1:3)
    isnothing(band) && throw(DomainError(x,"x is outside the working bands"))
    wall=band==1 ? -1 : 1
    v=_weld_height(target,x)
    v1,v2=if hasproperty(solution,:seam)
        wall == -1 ? solution.seam.left : solution.seam.right
    else
        map(i -> welded_wall_height(f,i,wall,f.m),(1,2))
    end
    if v<=v1 || wall== -1 || v>=v2
        i=v<=v1 ? 1 : 2
        t=_weld_tile_at_height(f,i,wall,v)
        b=f.blocks[i]
        z=_weighted_phase(b.c,b.T,b.seed,t+f.offset,wall)
        width=wall==1 ? b.widthfactor : 1.
        y=wall+f.epsilon*width*z
        return (;value=y,normalized=(band==3 ? f.third_weight : 1.)*width*z)
    end
    # The exposed part of the narrow rectangle's horizontal side belongs
    # to the real boundary and can lie inside the broad working interval.
    function connector_height(s)
        s==f.edge_end && return welded_wall_height(f,f.bi,1,f.m)
        z=_weld_remote(f,f.ai,s,f.h[f.ai])
        f.ai==1 ? z : f.modulus-z
    end
    s=_weld_bisect(connector_height,v,f.edge_end,1.,f.bi==1)
    a=f.blocks[f.ai]
    z=(1-a.delta)*real(elliptic_value(2s-1,a.T))
    y=(a.edge[1,1]*z+a.edge[1,2])/(a.edge[2,1]*z+a.edge[2,2])
    (;value=y,normalized=(y-1)*(band==3 ? f.third_weight : 1.)/f.epsilon)
end
