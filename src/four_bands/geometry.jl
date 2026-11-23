"""Schwarz–Christoffel differential and the measured heights of its four walls."""
struct Geometry
    e::Vector{Float64}
    p::Vector{Float64}
    H::Vector{Float64}
    roots::Vector{Float64}
    slit::Vector{Float64}
    left::Float64
    order::Int
    mom::Matrix{Float64}
end

struct FilterModel
    geom::Geometry
    endpoints::Vector{Float64}
    n::Int
    m1::Int
    m2::Int
    T::Float64
    trims::NTuple{2,Float64}
    sides::NTuple{2,Symbol}
    diagnostics::NamedTuple
end

struct GeometryError <: Exception
    message::String
end
Base.showerror(io::IO, e::GeometryError) = print(io, e.message)

const _gauss_rules = Dict{Int,Tuple{Vector{Float64},Vector{Float64}}}()
const _gauss_lock = ReentrantLock()

function _gauss(order::Int)
    order >= 2 || throw(ArgumentError("quadrature order must be at least two"))
    lock(_gauss_lock) do
        get!(_gauss_rules, order) do
            # Golub-Welsch for Legendre polynomials on [-1, 1].
            offdiag = [j / sqrt(4.0*j^2-1) for j in 1:order-1]
            eig = eigen(SymTridiagonal(zeros(order), offdiag))
            (eig.values, 2 .* vec(eig.vectors[1,:]).^2)
        end
    end
end

function _check_endpoints(e)
    length(e) == 8 && all(isfinite,e) && all(>(0),diff(e)) ||
        throw(ArgumentError("expected eight finite, strictly increasing endpoints"))
end

"""Moments of 1,x,x² against dx/|y| on interval i; optional angle limits."""
function moments(e, i, order=128, lo=0.0, hi=π)
    nodes, weights = _gauss(order)
    a, b = e[i], e[i+1]
    mid, half = (Float64(lo)+Float64(hi))/2, (Float64(hi)-Float64(lo))/2
    sums = zeros(3)
    correction = zeros(3)
    @inbounds for k in eachindex(nodes)
        angle = mid + half*nodes[k]
        x = a + (b-a)*cos(angle/2)^2
        product = 1.0
        for j in eachindex(e)
            if j != i && j != i+1
                product *= abs(x-e[j])
            end
        end
        v = half*weights[k]/sqrt(product)
        for j in 1:3
            increment = v - correction[j]
            next = sums[j]+increment
            correction[j] = (next-sums[j])-increment
            sums[j] = next
            v *= x
        end
    end
    sums
end

_angle(a,b,x) = 2atan(sqrt(b-x),sqrt(x-a))
_poly(p,x) = evalpoly(x,(p[1],p[2],p[3]))

function geometry(endpoints, order::Int=128)
    e = Float64.(collect(endpoints))
    _check_endpoints(e)
    M = Matrix{Float64}(undef,7,3)
    for i in 1:7
        M[i,:] = moments(e,i,order)
    end
    p = M[[2,4,6],:] \ [-2.0,0.0,0.0]
    d = muladd(-4p[1],p[3],p[2]^2)
    d > 0 && p[3] < 0 || throw(GeometryError("ill-conditioned SC quadratic"))
    q = -(p[2] + copysign(sqrt(d),p[2]))/2
    roots = sort([q/p[3],p[1]/q])
    slit = zeros(2)
    for (j,i) in enumerate((4,6))
        e[i] < roots[j] < e[i+1] || throw(GeometryError("SC root outside its slit gap"))
        int = dot(p,moments(e,i,order,_angle(e[i],e[i+1],roots[j]),π))
        slit[j] = 1 + (j == 1 ? int : -int)
    end
    H = (M[[3,5,7],:]*p) .* [-1.,1.,-1.]
    left = -dot(view(M,1,:),p)
    all(h->isfinite(h) && h>0,H) && isfinite(left) && left>0 ||
        throw(GeometryError("invalid SC wall heights"))
    Geometry(e,p,H,roots,slit,left,order,M)
end

function _wall_height(g, x, band)
    1 <= band <= 4 || throw(ArgumentError("band index must be in 1:4"))
    i = 2band-1
    a, b = g.e[i:i+1]
    a <= x <= b || throw(ArgumentError("point outside the selected band"))
    start = band == 1 ? g.left : sum(view(g.H,1:band-2))
    increment = dot(g.p,moments(g.e,i,g.order,_angle(a,b,x),π))
    start + (isodd(band) ? increment : -increment)
end

height_tiles(m::FilterModel,x::Real,band::Integer) = _wall_height(m.geom,x,band)/m.T

function _rule(f,a,b,order)
    nodes, weights = _gauss(order)
    mid, half = (a+b)/2, (b-a)/2
    half*sum(weights[i]*f(mid+half*nodes[i]) for i in eachindex(nodes))
end

function _quad(f,a=0.,b=1.,tol=1e-12,depth=18)
    low, high = _rule(f,a,b,16), _rule(f,a,b,32)
    abs(low-high) <= tol && return high
    depth > 0 || throw(GeometryError("adaptive SC quadrature did not converge"))
    mid = (a+b)/2
    _quad(f,a,mid,tol/2,depth-1) + _quad(f,mid,b,tol/2,depth-1)
end

function _upper_coordinate(g,z)
    base = g.e[3]
    function from_base(w,tol)
        delta = w-base
        _quad(0.,1.,tol) do t
            x = base + delta*t^2
            2sqrt(delta)*_poly(g.p,x)/prod(sqrt(x-g.e[j]) for j in eachindex(g.e) if j != 3)
        end
    end
    try
        return 1+from_base(z,1e-12)
    catch err
        err isa GeometryError || rethrow()
    end
    # Near another endpoint the direct segment can hug the real boundary.
    h = max(imag(z),abs(real(z)-base)/2,min(base-g.e[2],g.e[4]-base))
    a,b = complex(base,h),complex(real(z),h)
    function segment(a,b)
        a==b && return 0im
        _quad(0.,1.,1e-12/3) do t
            x = a+(b-a)*t
            (b-a)*_poly(g.p,x)/prod(sqrt(x-e) for e in g.e)
        end
    end
    1+from_base(a,1e-12/3)+segment(a,b)+segment(b,z)
end

function _outer_abscissa(g,x)
    right = x > last(g.e)
    start, direction = right ? (last(g.e),1) : (first(g.e),-1)
    others = right ? view(g.e,1:length(g.e)-1) : view(g.e,2:length(g.e))
    s = last(g.e)-first(g.e)
    t = sqrt(abs(x-start))/sqrt(s)
    upper = t > 1 ? 1-1/(1+t) : t/(1+t)
    # Compactifying each outer ray by x=start±s(t/(1-t))² cancels the
    # endpoint singularity and gives a bounded integrand at infinity.
    direction + _quad(0.,upper) do u
        v = 1-u
        d = direction*s*u^2
        xn = start*v^2+d
        p = if length(g.e)==6
            (right ? 1 : -1)*(g.p[1]*v^2+g.p[2]*xn)
        else
            g.p[1]*v^4+g.p[2]*xn*v^2+g.p[3]*xn^2
        end
        den = prod(sqrt(abs((start-b)*v^2+d)) for b in others)
        direction*2sqrt(s)*v*p/den
    end
end

"""Normalized SC coordinate. In the lower half-plane return its reflected value."""
function sc_coordinate(m::FilterModel,z::Number)
    isfinite(z) || throw(ArgumentError("SC evaluation requires a finite point"))
    imag(z)>0 && return _upper_coordinate(m.geom,ComplexF64(z))
    imag(z)<0 && return conj(_upper_coordinate(m.geom,conj(ComplexF64(z))))
    g, x = m.geom,Float64(real(z))
    for band in 1:4
        i = 2band-1
        if g.e[i] <= x <= g.e[i+1]
            return complex(band==1 ? -1. : 1.,_wall_height(g,x,band))
        end
    end
    for (i,level) in zip((2,4,6),(0,m.m1,m.m2))
        if g.e[i] < x < g.e[i+1]
            inc = dot(g.p,moments(g.e,i,2g.order,_angle(g.e[i],g.e[i+1],x),π))
            return complex((i==2 ? -1. : 1.)+(i==4 ? inc : -inc),level*m.T)
        end
    end
    complex(_outer_abscissa(g,x),m.n*m.T)
end
