# Genus-three periods and theta inversion on the marked real hyperelliptic curve.
# The cutoff controls a finite lattice sum in each reduced coordinate. The
# residual tolerances below are numerical diagnostics in Float64.
using LinearAlgebra

const BRANCH_ALPHA = [0 0 0; 1 0 0; 1 0 0; 0 1 0; 0 1 0; 0 0 1; 0 0 1; 0 0 0]
const BRANCH_BETA = [0 0 0; 0 0 0; 1 0 0; 1 0 0; 1 1 0; 1 1 0; 1 1 1; 1 1 1]
const RIEMANN_A = [1, 1, 1]
const RIEMANN_B = [1, 0, 1]

theta_char(indices...) = (
    mod.(reduce(+, (BRANCH_ALPHA[i, :] for i in indices); init=RIEMANN_A), 2),
    mod.(reduce(+, (BRANCH_BETA[i, :] for i in indices); init=RIEMANN_B), 2)
)

struct ThetaCurve
    e::Vector{Float64}
    omega::Matrix{Float64}
    coefficients::Matrix{Float64}
    branch_aj::Matrix{ComplexF64}
    cutoff::Int
end

"""Build periods from the seven rows of actual gap and band moments."""
function ThetaCurve(g; cutoff=nothing)
    m = g.mom
    A = 2 .* m[[2, 4, 6], :]
    A[[1, 3], :] .*= -1
    B = -2 .* vcat(m[1:1, :], m[1:1, :] - m[3:3, :], m[1:1, :] - m[3:3, :] + m[5:5, :])
    raw = B / A
    maximum(abs, raw - transpose(raw)) <= 1e-8 ||
        throw(ArgumentError("SC period quadrature does not give a symmetric matrix"))
    Ω = (raw + transpose(raw)) / 2
    λ = eigmin(Symmetric(Ω))
    λ > 0 || throw(ArgumentError("SC period matrix is not positive definite"))
    # A Gaussian decay estimate chooses an initial cutoff. Condition numbers,
    # period quadrature and floating-point errors can require extra precision.
    N = isnothing(cutoff) ? max(7, ceil(Int, 1 + sqrt(40/(π*λ)))) : cutoff
    N isa Integer && !(N isa Bool) && 1 <= N <= 64 ||
        throw(ArgumentError("theta cutoff must be between 1 and 64; check period conditioning"))
    ThetaCurve(copy(g.e), Ω, inv(A), (BRANCH_BETA - im*BRANCH_ALPHA*Ω)/2, N)
end

"""Characteristic theta and its complex gradient at a reduced argument."""
function theta_value_gradient(c::ThetaCurve, z, characteristic; scaled=false)
    Ω = c.omega
    ai, bi = characteristic
    a, b = mod.(ai, 2), mod.(bi, 2)
    charsign = isodd(dot(a, (bi .- b) .÷ 2)) ? -1 : 1
    shifts = Ω \ imag.(z)
    q = round.(shifts)
    p = round.(real.(z))
    r = real.(z) - p
    phase = π * (dot(a, p) - dot(b, q) - 2*dot(q, r))
    shift = shifts - q
    v = r + b/2
    # Completing the square gives the scaled Gaussian directly, even on
    # translated arguments whose unscaled multiplier would overflow.
    value = 0.0im
    g1 = g2 = g3 = 0.0im
    N = c.cutoff
    @inbounds for k3 in -N:N, k2 in -N:N, k1 in -N:N
        x, y, t = k1+a[1]/2, k2+a[2]/2, k3+a[3]/2
        X, Y, Z = x+shift[1], y+shift[2], t+shift[3]
        quadratic = Ω[1,1]*X^2 + Ω[2,2]*Y^2 + Ω[3,3]*Z^2 +
                    2*(Ω[1,2]*X*Y + Ω[1,3]*X*Z + Ω[2,3]*Y*Z)
        term = exp(-π*quadratic) * cis(2π*(x*v[1] + y*v[2] + t*v[3]))
        value += term
        g1 += x*term
        g2 += y*term
        g3 += t*term
    end
    factor = charsign * cis(phase)
    if !scaled
        factor *= exp(π*dot(imag.(z), shifts))
    end
    isfinite(factor) || throw(DomainError(z, "theta multiplier exceeds Float64 range"))
    return factor*value, (2π*im*factor) .* ([g1, g2, g3] - q*value)
end

riemann_constant_residual(c::ThetaCurve) =
    abs(first(theta_value_gradient(c, zeros(3), (RIEMANN_A, RIEMANN_B))))
theta_null_residual(c::ThetaCurve) = riemann_constant_residual(c)

function projection_chart(c::ThetaCurve, seed)
    # Mobius normalization at branch points 1, 8 and 2. Choose an auxiliary
    # pair away from common zeros, and keep this pair fixed in one Newton solve.
    j = c.branch_aj[2, :]
    best = 0.0
    chart = nothing
    for a in 3:6, b in a+1:7
        cs, cl = theta_char(1, a, b), theta_char(8, a, b)
        numerator = first(theta_value_gradient(c, j, cs))
        denominator = first(theta_value_gradient(c, j, cl))
        abs(numerator*denominator) < 1e-24 && continue
        d = abs(first(theta_value_gradient(c, seed, cl)))
        if d > best && d > 1e-14
            best = d
            chart = (cs, cl, (denominator/numerator)^2)
        end
    end
    isnothing(chart) && throw(DomainError(seed, "all projection charts are singular"))
    return chart
end

function theta_system(c, z, target, chart)
    cs, cl, constant = chart
    r = Vector{ComplexF64}(undef, 5)
    J = Matrix{ComplexF64}(undef, 5, 3)
    for i in 1:4
        r[i], J[i, :] = theta_value_gradient(c, z, theta_char(i))
    end
    top, gt = theta_value_gradient(c, z, cs)
    bottom, gb = theta_value_gradient(c, z, cl)
    abs(bottom) > 1e-14 || throw(DomainError(z, "projection denominator vanishes"))
    ratio = top/bottom
    s = max(1, abs(target))
    r[5] = (constant*ratio^2 - target)/s
    J[5, :] = (2*constant*ratio/(bottom*s)) .* (gt - ratio .* gb)
    return r, J
end

function theta_correct(c, x, seed)
    e = c.e
    target = (x-e[1])/(x-e[8]) * (e[2]-e[8])/(e[2]-e[1])
    chart = projection_chart(c, seed)
    z = copy(seed)
    r, J = theta_system(c, z, target, chart)
    for _ in 1:32
        maximum(abs, r) < 2e-12 && return z
        step = J \ r
        accepted = false
        for j in 0:12
            next = z - exp2(-j)*step
            trial = try
                theta_system(c, next, target, chart)
            catch err
                err isa DomainError || rethrow()
                continue
            end
            next_r, next_J = trial
            if norm(next_r) < norm(r)
                z, r, J = next, next_r, next_J
                accepted = true
                break
            end
        end
        accepted || break
    end
    maximum(abs, r) <= 5e-10 && return z
    throw(DomainError(x, "theta inversion residual exceeds tolerance"))
end

"""Evaluate the marked SC coordinate by continuation of four theta loci and the projection."""
function theta_coordinate(c::ThetaCurve, point::Number)
    x = ComplexF64(point)
    isfinite(x) && imag(x) > 0 || throw(ArgumentError("theta continuation needs Im(x)>0"))
    base = c.e[3]
    delta = x-base
    gap = min(base-c.e[2], c.e[4]-base)
    t = min(0.01*gap/abs(delta), 1.0)
    at = base+t*delta
    root = prod(sqrt(complex(base-c.e[i], 0.0)) for i in 1:8 if i != 3)
    transform = transpose(c.coefficients)
    seed = c.branch_aj[3, :] + (2*sqrt(at-base)/root) .* (transform * [1, base, base^2])
    z = theta_correct(c, at, seed)
    dt = t
    while t < 1
        next_t = min(t+dt, 1.0)
        next_x = base + next_t*delta
        dz = transform * [1, at, at^2] / prod(sqrt(at-e) for e in c.e)
        next_z = try
            theta_correct(c, next_x, z + (next_x-at)*dz)
        catch err
            err isa DomainError || rethrow()
            dt /= 2
            (dt < 1e-9 || t+dt == t) && rethrow()
            continue
        end
        z, t, at = next_z, next_t, next_x
        dt = min(1.7*dt, 1-t)
    end
    return 1 + 4*(z[1]-c.branch_aj[3, 1])
end

"""An exact LDL test for a positive lower eigenvalue bound of the stored matrix."""
function theta_spectral_bound(c::ThetaCurve)
    λ = eigmin(Symmetric(c.omega))/2
    for _ in 1:64
        λ > 0 && isfinite(λ) || break
        A = Rational{BigInt}.(c.omega) - Rational{BigInt}(λ)*I
        valid = true
        for k in 1:3
            pivot = A[k,k]
            if pivot <= 0
                valid = false
                break
            end
            for i in k+1:3, j in i:3
                A[i,j] -= A[i,k]*A[j,k]/pivot
                A[j,i] = A[i,j]
            end
        end
        valid && return λ
        λ /= 2
    end
    throw(DomainError(c.omega, "could not certify a positive spectral lower bound"))
end

# Bound by a union over the three coordinates outside [-N,N] and two
# complete one-dimensional Gaussian sums. Here d >= max |shift[i]|.
function _theta_logtail(λ, d, N, logfactor)
    return logfactor + log(3) + 2log1p(inv(sqrt(λ))) - 0.5log(λ) - π*λ*(N-d)^2
end

"""Choose a cutoff from Gaussian tail bounds for the scaled theta and gradient.

The bounds apply to series truncation on the stored matrix. The LDL test
uses exact rational arithmetic on its Float64 entries. Quadrature, other
roundoff, and nonlinear inversion are outside this error estimate.
"""
function theta_controlled(c::ThetaCurve, z, characteristic; tolerance=1e-12, scaled=true)
    length(z) == 3 && all(isfinite,z) || throw(ArgumentError("expected three finite coordinates"))
    length(characteristic) == 2 && all(v->length(v)==3 && all(isinteger,v), characteristic) ||
        throw(ArgumentError("expected two length-three integer characteristics"))
    isfinite(tolerance) && tolerance > 0 || throw(ArgumentError("expected a positive tolerance"))
    λ = theta_spectral_bound(c)
    shifts = c.omega \ imag.(z)
    d = maximum(abs, shifts - round.(shifts) + mod.(first(characteristic),2)/2)
    logfactor = scaled ? 0.0 : π*dot(imag.(z), shifts)
    # In a derivative term |l_i-q_i| <= ||l+shift|| + |shifts_i|. Absorb
    # its first summand using r*exp(-πλr²/2) <= 1/sqrt(πλe).
    gradfactor = 2π * (inv(sqrt(π*λ*exp(1))) + maximum(abs, shifts))
    for N in max(1, floor(Int,d)+1):64
        value_tail = _theta_logtail(λ,d,N,logfactor)
        grad_tail = _theta_logtail(λ/2,d,N,logfactor) + log(gradfactor)
        if max(value_tail,grad_tail) <= log(tolerance)
            localcurve = ThetaCurve(c.e,c.omega,c.coefficients,c.branch_aj,N)
            value, gradient = theta_value_gradient(localcurve,z,characteristic; scaled=scaled)
            return (; value, gradient, cutoff=N, lambda_lower=λ,
                      log_value_tail=value_tail, log_gradient_tail=grad_tail)
        end
    end
    throw(DomainError(z, "theta truncation tolerance needs a cutoff above 64"))
end
