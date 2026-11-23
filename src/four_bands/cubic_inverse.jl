"""Degree-three all-endpoint alternant for targets −,+,+,+ (canonical source genus five)."""
struct CubicFourBands
    endpoints::Vector{Float64}
    n::Int
    T::Float64
    mu::Float64
    log_mu::Float64
    branches::NTuple{8,BigFloat}
    deviation::BigFloat
    lower::BigFloat
    upper::BigFloat
end

function _cubic_agm(x)
    a,b = one(x),x
    for _ in 1:64
        abs(a-b)<=4eps(BigFloat)*a && return (a+b)/2
        a,b = (a+b)/2,sqrt(a*b)
    end
    throw(GeometryError("elliptic arithmetic-geometric mean failed to converge"))
end

"""
    cubic_fourbands(endpoints)

Construct the explicit real rational cubic alternating at all eight prescribed
endpoints (the first interval has target −1). Its companion is (1−μ²)/r.
Eight finite, strictly increasing Float64 endpoints are required. Products and
cross ratios use extended precision before reporting Float64 values.
"""
function cubic_fourbands(endpoints)
    e = Float64.(collect(endpoints))
    _check_endpoints(e)
    t = Tuple(BigFloat.(e))
    x,y = t[1],t[2]
    # For each pair the factor of ρ=H(y)/H(x) is 1−δ. Its direct
    # cross-ratio form resolves ρ≈0; log1p(−δ) handles ρ≈1.
    logρ = zero(x)
    for j in (3,5,7)
        a,b = t[j],t[j+1]
        den = (b-y)*(a-x)
        δ = (y-x)*(b-a)/den
        logρ += δ<.5 ? log1p(-δ) : log((a-y)*(b-x)/den)
    end
    μ = sqrt(-expm1(logρ))
    U = 1+μ
    L = exp(logρ)/U
    # K(k')/K(k) with k=L/U and k'=2√μ/U. Independent complements
    # avoid subtracting two nearly equal numbers at small or large μ.
    T = _cubic_agm(2sqrt(μ)/U)/_cubic_agm(L/U)
    CubicFourBands(e,3,Float64(T),Float64(μ),Float64(log(μ)),t,μ,L,U)
end

# A and B are homogeneous numerator and denominator for H(z)/H(E₁).
# D=B−A is evaluated by telescoping the products of linear factors.
function _cubic_products(f,z)
    t = f.branches
    x = t[1]
    A = one(z)
    B = one(z)
    D = zero(z)
    for j in (3,5,7)
        a,b = t[j],t[j+1]
        af = (z-a)/(x-a)
        bf = (z-b)/(x-b)
        difference = (z-x)*(b-a)/((x-a)*(x-b))
        D = D*bf+A*difference
        A *= af
        B *= bf
    end
    A,B,D
end

function _cubic_pair(f,z,(A,B,D)=_cubic_products(f,z))
    μ,L,U = f.deviation,f.lower,f.upper
    # When μ is small D avoids cancellation at the stopband. When μ≈1,
    # L is retained independently to resolve numerator zeros near a passband.
    num,den = μ<.25 ? (D-μ*B,D+μ*B) : (L*B-A,U*B-A)
    iszero(num) && throw(DomainError(z,"reciprocal cubic pole"))
    iszero(den) && throw(DomainError(z,"direct cubic pole"))
    U*num/den,L*den/num
end

function wall_values(f::CubicFourBands,x::Real,b::Integer)
    1<=b<=4 || throw(ArgumentError("expected a band in 1:4"))
    lo,hi = f.endpoints[2b-1:2b]
    lo<=x<=hi || throw(ArgumentError("point outside selected cubic band"))
    target = b==1 ? -1 : 1
    μ,L,U = f.deviation,f.lower,f.upper
    if x==lo || x==hi
        ripple = x==lo ? -1. : 1.
        low = target<0 ? -U : L
        high = target<0 ? -L : U
        a,c = x==lo ? (low,high) : (high,low)
        return (;direct=Float64(a),reciprocal=Float64(c),ripple,
                 reciprocal_ripple=-ripple,target,method=:rational)
    end
    z=BigFloat(x)
    A,B,D = _cubic_products(f,z)
    a,c = _cubic_pair(f,z,(A,B,D))
    if target>0
        ripple = -(U*B+A)/(U*B-A)
        reciprocal_ripple = (L*B+A)/(L*B-A)
    elseif μ<.25
        # 0≤s≤1 even when μ² underflows Float64.
        s = D/(μ^2*B)
        ripple = ((2+μ)*s-1)/(1+μ*s)
        reciprocal_ripple = (1-(2-μ)*s)/(1-μ*s)
    else
        ripple = (a+1)/μ
        reciprocal_ripple = (c+1)/μ
    end
    (;direct=Float64(a),reciprocal=Float64(c),ripple=Float64(ripple),
      reciprocal_ripple=Float64(reciprocal_ripple),target,method=:rational)
end

function evaluate(f::CubicFourBands,z::Number;method=:rational)
    method===:rational || throw(ArgumentError("the explicit cubic uses method=:rational"))
    isfinite(z) || throw(ArgumentError("expected a finite point"))
    if isreal(z)
        x=real(z)
        for b in 1:4
            f.endpoints[2b-1]<=x<=f.endpoints[2b] && return wall_values(f,x,b)
        end
        a,c=_cubic_pair(f,BigFloat(x))
        return (;direct=Float64(a),reciprocal=Float64(c),method)
    end
    a,c=_cubic_pair(f,Complex{BigFloat}(z))
    (;direct=ComplexF64(a),reciprocal=ComplexF64(c),method)
end

_plot_counts(::CubicFourBands) = (1,1,1,1)
_plot_description(::CubicFourBands) = "Дробь степени 3 по восьми концам · род 5"
_plot_amplitude(f::CubicFourBands) = (;mu=f.mu,log_mu=f.log_mu)
