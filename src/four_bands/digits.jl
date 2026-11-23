# The entries α, β are binary columns in the convention of (14) and (19).
# A digit denotes a marked branch point; the Riemann constant is itself a word.
function _branch_indices(word::AbstractString)
    s = strip(word)
    if startswith(s, "[") && endswith(s, "]")
        s = chop(s; head=1, tail=1)
    end
    all(c -> '1' <= c <= '9', s) || throw(ArgumentError("expected branch digits such as 35 or [35]"))
    return [Int(c)-Int('0') for c in s]
end

_branch_indices(word::Integer) = _branch_indices(string(word))
_branch_indices(word::Union{Tuple,AbstractVector}) = word

"""
    branch_characteristic(g, word)

Return the binary columns (α, β) of the characteristic [word], the sum
modulo two of the Abel–Jacobi images of the listed branch points, based at
point 1. Examples include `branch_characteristic(2,"35")` for the Riemann
constant of the article and `branch_characteristic(3,"357")` for its
counterpart in the project's eight-point marking. Use explicit index lists
for multi-digit branch labels, e.g. `(3,5,7,9,11)` in genus five.
"""
function branch_characteristic(g::Integer, word)
    1 <= g <= 5 || throw(ArgumentError("branch digit notation supports genera 1 to 5"))
    a, b = zeros(Int,g), zeros(Int,g)
    for s in _branch_indices(word)
        s isa Integer && !(s isa Bool) && 1 <= s <= 2g+2 ||
            throw(ArgumentError("branch index must be in 1:$(2g+2)"))
        j = s ÷ 2
        if 1 <= j <= g
            a[j] = 1-a[j]
        end
        for i in 1:((s-1) ÷ 2)
            b[i] = 1-b[i]
        end
    end
    return a, b
end

# The real part of Π may be nonzero. Complete the Gaussian square and
# re-index the sum to keep its peak inside the truncation box.
function _digit_reduction(P,z,a,b)
    X, Y = real(P), imag(P)
    v = imag(z)
    t = Y \ v
    q = round.(t)
    sh = t-q
    h = real(z) - X*q
    p = round.(h)
    h = h-p+b/2
    factor = exp(π*dot(v,t))*cispi(dot(q,X*q)-2dot(q,real(z)+b/2)+dot(a,p))
    isfinite(factor) || throw(DomainError(z, "unscaled theta exceeds Float64 range"))
    return X,Y,sh,h,factor
end

function _digit_sum(P, z, a, b, N, ::Val{G}) where G
    X,Y,sh,h,factor=_digit_reduction(P,z,a,b)
    result = 0.0im
    ranges = ntuple(i -> -N:N, Val(G))
    for k in CartesianIndices(ranges)
        x = ntuple(i -> k[i]+a[i]/2, Val(G))
        d = ntuple(i -> x[i]+sh[i], Val(G))
        decay, phase = 0., 0.
        @inbounds for i in 1:G
            decay += Y[i,i]*d[i]^2
            phase += X[i,i]*x[i]^2 + 2x[i]*h[i]
            for j in 1:i-1
                decay += 2Y[i,j]*d[i]*d[j]
                phase += 2X[i,j]*x[i]*x[j]
            end
        end
        result += exp(-π*decay)*cispi(phase)
    end
    value = factor*result
    isfinite(value) || throw(DomainError(z,"theta overflows Float64"))
    return value
end

# For large g enumerate a Gaussian ellipsoid instead of its enclosing cube.
# Outside dᵀYd ≤ R² the Gaussian is bounded by exp(−πR²/2) times a
# second Gaussian with matrix Y/2, whose entire sum is bounded via λₘᵢₙ.
function _digit_ellipsoid(P,z,a,b,λ)
    X,Y,sh,h,factor=_digit_reduction(P,z,a,b)
    U=cholesky(Symmetric(Y)).U
    g=length(a)
    radius2=2*(40+g*log1p(sqrt(2/λ)))/π
    x=zeros(g)
    d=zeros(g)
    count=Ref(0)
    function accumulate(i,budget)
        if i==0
            count[]+=1
            count[]<=4_000_000 || throw(ArgumentError("theta ellipsoid exceeds the numerical work limit"))
            phase=0.
            @inbounds for k in 1:g
                phase+=X[k,k]*x[k]^2+2x[k]*h[k]
                for j in 1:k-1
                    phase+=2X[k,j]*x[k]*x[j]
                end
            end
            return exp(-π*(radius2-budget))*cispi(phase)
        end
        cross=0.
        @inbounds for j in i+1:g
            cross+=U[i,j]*d[j]
        end
        center=-(cross/U[i,i]+a[i]/2+sh[i])
        width=sqrt(budget)/U[i,i]
        isfinite(width) && width<=4_000_000 || throw(ArgumentError("theta ellipsoid needs better conditioning"))
        result=0.0im
        for k in ceil(Int,center-width):floor(Int,center+width)
            x[i]=k+a[i]/2
            d[i]=x[i]+sh[i]
            left=budget-(U[i,i]*d[i]+cross)^2
            left>=0 && (result+=accumulate(i-1,left))
        end
        return result
    end
    value=factor*accumulate(g,radius2)
    isfinite(value) || throw(DomainError(z,"theta overflows Float64"))
    return value
end

"""
    theta_digits(Π, word, z=zeros(size(Π,1)); cutoff=nothing)
    theta_digits(curve::ThetaCurve, word, z=zeros(3))

Evaluate the Fourier series (19) with the reduced characteristic [word],
in Float64. Π denotes the complex Riemann matrix, e.g. `im * Ω`.
The matrix overload evaluates theta on any symmetric matrix of the
supported size with positive definite imaginary part. Branch-point
interpretation additionally requires the corresponding marked Jacobian.
"""
function theta_digits(Pi::AbstractMatrix, word, z=zeros(size(Pi,1)); cutoff=nothing)
    g = size(Pi,1)
    a,b = branch_characteristic(g,word)
    P = Matrix{ComplexF64}(Pi)
    size(P,2)==g && all(isfinite,P) && issymmetric(P) ||
        throw(ArgumentError("Π must be finite, square and symmetric"))
    length(z)==g && all(isfinite,z) || throw(ArgumentError("expected $g finite Jacobi coordinates"))
    λ = eigmin(Symmetric(imag(P)))
    isfinite(λ) && λ>0 || throw(ArgumentError("Im Π must be positive definite; use Π = im * Ω"))
    # A reduced Gaussian has shifts of modulus at most one in each
    # coordinate. Budget for all coordinate tails at once (scaled series).
    target = 40 + log(2g) + g*log1p(inv(sqrt(λ)))
    radius = 1+sqrt(target/(π*λ))
    if isnothing(cutoff)
        g>=4 && return _digit_ellipsoid(P,ComplexF64.(collect(z)),a,b,λ)
        radius <= 64 || throw(ArgumentError("theta truncation needs more than 64 shells"))
        N = ceil(Int,radius)
    else
        cutoff isa Integer && !(cutoff isa Bool) && 1 <= cutoff <= 64 ||
            throw(ArgumentError("cutoff must be an integer in 1:64"))
        N = Int(cutoff)
        (2N+1)^g<=12_000_000 || throw(ArgumentError("explicit theta cutoff exceeds the numerical work limit; omit cutoff for ellipsoid summation"))
    end
    return _digit_sum(P, ComplexF64.(collect(z)), a, b, N, Val(g))
end

function theta_digits(c::ThetaCurve, word, z=zeros(3))
    length(z)==3 && all(isfinite,z) || throw(ArgumentError("expected three finite Jacobi coordinates"))
    return first(theta_value_gradient(c,ComplexF64.(collect(z)),branch_characteristic(3,word)))
end

struct _ThetaNotation end
struct _DigitTheta{W}
    word::W
end

"""Article notation: `θ[35](u, Π)`, `θ[35](Π)` or equivalently `θ[3,5](u, Π)`."""
const θ = _ThetaNotation()
Base.getindex(::_ThetaNotation) = _DigitTheta(())
Base.getindex(::_ThetaNotation,word) = _DigitTheta(word)
Base.getindex(::_ThetaNotation,i::Integer,j::Integer,ks::Integer...) = _DigitTheta((i,j,ks...))
(t::_DigitTheta)(z, P) = theta_digits(P,t.word,z)
(t::_DigitTheta)(P) = theta_digits(P,t.word)

"""
    paper_period_matrix(τ, n, m, Ω22)

Genus-two Π of the article (17), in its marked real chamber (12).
Its lower-right chamber bound is Ω22 > Ω12; §2.7 also prints a stronger
inequality immediately before describing the full interval (Ω12,∞).
"""
function paper_period_matrix(τ::Number, n::Integer, m::Integer, Ω22::Real)
    isfinite(τ) && iszero(real(τ)) && imag(τ)>0 ||
        throw(ArgumentError("τ must be purely imaginary with Im τ > 0"))
    1 <= m < n || throw(ArgumentError("expected 1 ≤ m < n"))
    T = Float64(imag(τ))
    o11, o12 = n*T/2, (n-m)*T/2
    isfinite(o11) && 0 < o12 < min(o11,Ω22) && isfinite(Ω22) ||
        throw(ArgumentError("expected a finite period matrix in the chamber of (12)"))
    return im * [o11 o12; o12 Float64(Ω22)]
end

"""
    paper_branchpoints(Π)

The six ordered real branch values (0,1,x₃,x₄,x₅,∞) of the paper's
genus-two projection. Computes (23) with digit characteristics from §2.3.
"""
function paper_branchpoints(Π::AbstractMatrix)
    size(Π)==(2,2) && all(iszero,real(Π)) || throw(ArgumentError("expected the purely imaginary 2×2 period matrix (11)"))
    Ω = imag(Π)
    0 < Ω[1,2] < min(Ω[1,1],Ω[2,2]) || throw(ArgumentError("period matrix outside chamber (12)"))
    # In (23), read both *columns* of each printed characteristic. For
    # example the rows (00;01) mean α=(0,0), β=(0,1), i.e. [45].
    t0,t01,t10,t11 = real(θ[1](Π)),real(θ[45](Π)),real(θ[23](Π)),real(θ[6](Π))
    taa,tab = real(θ[34](Π)),real(θ[25](Π))
    x3 = (t0*t01/(t10*t11))^2
    x4 = (t01*taa/(t10*tab))^2
    x5 = (t0*taa/(tab*t11))^2
    all(isfinite,(x3,x4,x5)) && 1 < x3 < x4 < x5 ||
        throw(DomainError(Π,"the branch points cannot be resolved with Float64"))
    return [0.,1.,x3,x4,x5,Inf]
end
