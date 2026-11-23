"""Elliptic modulus, complete integral and normalizations for τ=iT."""
function elliptic_constants(T::Real)
    T = Float64(T)
    isfinite(T) && T > 0 || throw(ArgumentError("T must be finite and positive"))
    if T < 0.5
        # Imaginary transformation. Keep θ₃−θ₄ as its own positive series.
        logq = -π/T
        isfinite(logq) || throw(DomainError(T, "log deviation is outside Float64"))
        odd = 0.0
        θ3c = θ4c = 1.0
        for j in 1:100
            term = exp(logq*j^2)
            θ3c += 2term
            θ4c += (isodd(j) ? -2 : 2) * term
            isodd(j) && (odd += term)
            term < eps(Float64)/4 && break
        end
        θ2, θ3 = θ4c/sqrt(T), θ3c/sqrt(T)
        k = (θ4c/θ3c)^2
        # Factor q out of μ so the logarithm stays accurate through underflow.
        mu = 4odd*(θ3c+θ4c)/(θ3c^2+θ4c^2)
        log_mu = logq < -36 ? log(4.0)+logq : log(mu)
    else
        θ3, θ2 = 1.0, 0.0
        for j in 0:100
            term = 2exp(-π*T*(j+0.5)^2)
            θ2 += term
            j > 0 && (θ3 += 2exp(-π*T*j^2))
            j > 2 && term < eps(Float64)*θ2 && break
        end
        k = (θ2/θ3)^2
        k > 0 || throw(DomainError(T, "elliptic modulus underflows Float64"))
        mu = (1-k)/(1+k)
        log_mu = log(mu)
    end
    K = (π/2)*θ3^2
    isfinite(K) || throw(DomainError(T, "elliptic integral is outside Float64"))
    return (; k, K, mu, log_mu, scale=2k/(1+k), theta2=θ2, theta3=θ3)
end

# In this cell |Re u|≤1, |Im u|≤T/2. Combine Gaussian damping with
# the imaginary argument before exponentiating on long rectangles.
function _cell_sn(u, T, c)
    numerator = 0.0im
    denominator = 0.0im
    if T < 0.5
        w = u/(2im*T)
        base = π*abs(imag(w))
        for j in 0:100
            a = (2j+1)*π*w
            damping = -π*j*(j+1)/T-base
            pos = exp(complex(damping-imag(a), real(a)))
            neg = exp(complex(damping+imag(a), -real(a)))
            sine, cosine = (pos-neg)/(2im), (pos+neg)/2
            numerator += (isodd(j) ? -sine : sine)
            denominator += cosine
            j > 1 && max(abs(sine), abs(cosine)) < eps(Float64)/4 && break
        end
        return im*(c.theta3/c.theta2)*numerator/denominator
    end
    denominator = 1.0 + 0.0im
    for j in 0:100
        # Divide θ₁ by exp(−πT/4) before summing; remove the same
        # factor from the theta constant θ₂ in its prefactor below.
        a = (j+0.5)*π*u
        damping = -π*T*j*(j+1)
        # For large T the reduced sine may grow, but k>0 keeps it finite.
        s = (exp(complex(damping-imag(a), real(a))) -
             exp(complex(damping+imag(a), -real(a))))/im
        numerator += isodd(j) ? -s : s
        if j > 0
            b = j*π*u
            d = -π*T*j^2
            cosine = exp(complex(d-imag(b), real(b))) +
                     exp(complex(d+imag(b), -real(b)))
            denominator += isodd(j) ? -cosine : cosine
            j > 2 && max(abs(s), abs(cosine)) < eps(Float64)/4 && break
        end
    end
    return (c.theta3/(c.theta2/exp(-π*T/4)))*numerator/denominator
end

# Real/imaginary translations give a sign and may swap the two companions.
function _elliptic_cell(u, T, c)
    isfinite(u) || throw(ArgumentError("finite SC coordinate required"))
    x = rem(real(u), 4.0, RoundNearest)
    y = rem(imag(u), 2T, RoundNearest)
    shiftx, shifty = abs(x)>1, abs(y)>T/2
    shiftx && (x -= copysign(2.0, x))
    shifty && (y -= copysign(T, y))
    raw = _cell_sn(complex(x,y), T, c)
    return (shiftx ? -raw : raw), shifty
end

"""Jacobi sn(Ku | k²), where T determines the modulus k."""
function elliptic_value(u, T::Real)
    c = elliptic_constants(T)
    raw, swapped = _elliptic_cell(u, T, c)
    swapped && iszero(raw) && throw(DomainError(u, "elliptic pole"))
    return swapped ? inv(c.k*raw) : raw
end

function _pair(u, T, c)
    raw, swapped = _elliptic_cell(u, T, c)
    iszero(raw) && throw(DomainError(u, swapped ? "direct pole" : "reciprocal pole"))
    a, b = c.scale*raw, (2/(1+c.k))/raw
    isfinite(a) && isfinite(b) || throw(DomainError(u, "unresolved pole or overflow"))
    return swapped ? (b,a) : (a,b)
end

function _wall_ripple(phase, T, c, target)
    if T >= 0.5
        direct, _ = _pair(1+im*phase*T, T, c)
        return target*(real(direct)-1)/c.mu
    end
    # Complementary sn, with q^(1/4) removed from θ₁ and θ₂.
    q = exp(-π/T)
    θ1, θ2 = 0.0, 0.0
    θ3, θ4 = 1.0, 1.0
    for j in 0:100
        term = q^(j*(j+1))
        s = term*sinpi((2j+1)*phase/2)
        θ1 += isodd(j) ? -s : s
        θ2 += term
        if j > 0
            v = q^(j^2)
            θ3 += 2v
            θ4 += (isodd(j) ? -2 : 2)*v*cospi(j*phase)
        end
        j > 1 && term < eps(Float64)/4 && break
    end
    sn = θ3*θ1/(θ2*θ4)
    dn = sqrt(1 - 16q*(θ2/θ3)^4*sn^2)
    return target*(2c.k*(1+c.k)*sn^2/(dn*(1+dn))-1)
end

"""Signed error in units of μ on a vertical wall, including underflow of μ."""
function boundary_ripple(height::Real, T::Real, target::Integer=1)
    abs(target)==1 || throw(ArgumentError("target must be ±1"))
    c = elliptic_constants(T)
    phase = rem(height/T, 2.0, RoundNearest)
    return _wall_ripple(phase,T,c,target)
end

"""Evaluate both fractions and their μ-scaled errors on maximal band 1…4."""
function wall_values(model, x::Real, band::Integer)
    1 <= band <= 4 || throw(ArgumentError("band must be 1, 2, 3 or 4"))
    _wall_values_at(height_tiles(model,x,band),model.T,band==1 ? -1 : 1)
end

function _wall_values_at(tiles,T,target)
    c = elliptic_constants(T)
    phase = rem(tiles, 2.0, RoundNearest)
    ripple = _wall_ripple(phase, T, c, target)
    reciprocal_ripple = _wall_ripple(phase+1, T, c, target)
    if c.mu < 0.5
        direct = target + c.mu*ripple
        reciprocal = target + c.mu*reciprocal_ripple
    else
        a,b = _pair(target + im*phase*T, T, c)
        direct,reciprocal = real(a),real(b)
    end
    return (; direct, reciprocal, ripple, reciprocal_ripple, target,
              height_tiles=tiles, method=:wall)
end

"""Evaluate both fractions, using :sc or :theta in either open half-plane."""
function evaluate(model, z::Number; method=:sc, curve=nothing)
    method in (:sc,:theta) || throw(ArgumentError("method must be :sc or :theta"))
    isfinite(z) || throw(ArgumentError("evaluation point must be finite"))
    e = model.geom.e
    T = model.T
    if isreal(z)
        x = real(z)
        for band in 1:4
            e[2band-1] <= x <= e[2band] && return wall_values(model,x,band)
        end
        u = sc_coordinate(model, x)
        # A horizontal SC boundary sits on an integer height. Evaluate at
        # height zero, then use the exact integer parity for the swap.
        level = x < e[1] || x > e[8] ? model.n :
            x < e[3] ? 0 : x < e[5] ? model.m1 : model.m2
        abs(real(u)) > 2e-12 || throw(DomainError(x,
            isodd(level) ? "direct pole or unresolved SC abscissa" :
                           "reciprocal pole or unresolved SC abscissa"))
        a,b = _pair(real(u), T, elliptic_constants(T))
        direct, reciprocal = isodd(level) ? (real(b),real(a)) : (real(a),real(b))
        return (; direct, reciprocal, height_tiles=level, sc_abscissa=real(u),
                  method=:horizontal)
    end
    upper = complex(real(z), abs(imag(z)))
    u = if method == :theta
        localcurve = isnothing(curve) ? ThetaCurve(model.geom) : curve
        theta_coordinate(localcurve, upper)
    else
        sc_coordinate(model, upper)
    end
    direct, reciprocal = _pair(u, T, elliptic_constants(T))
    if imag(z) < 0
        direct, reciprocal = conj(direct), conj(reciprocal)
    end
    return (; direct, reciprocal, method)
end
