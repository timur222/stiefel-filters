using HCubature
using QuadGK

const BIG_CONST = 1e5

const tol = 1e-6

function mesure_bounds(bound=BIG_CONST)
    return (
        -bound * 1.08,
        bound,
        -bound * 1.07,
        bound,
        -bound * 1.06,
        bound,
    )
end

function smal_mesure(y1, y2, y3, x)
    numerator = (y1 - y3) * (y2 - y3) * (y2 - y1)
    denominator = (x * (y2 - y3) + (y2 - y1))^2 + (x * y1 * (y2 - y3) + y3 * (y2 - y1))^2

    if iszero(denominator)
        return zero(promote_type(typeof(numerator), typeof(denominator)))
    end

    return numerator / denominator
end

function mesure(y1, y2, y3, x1, x2, x3)
    a1 = smal_mesure(y1, y2, y3, x1) / (1 + y1^2)
    a2 = smal_mesure(y1, y2, y3, x2) / (1 + y2^2)
    a3 = smal_mesure(y1, y2, y3, x3) / (1 + y3^2)
    return a1 * a2 * a3
end

function mesure_integral_y3(y1, y2, x1, x2, x3, y3_left, y3_right, tol)
    return quadgk(y3 -> mesure(y1, y2, y3, x1, x2, x3), y3_left, y3_right; rtol=tol)
end

function mesure_integral_y2(y1, x1, x2, x3, y2_left, y2_right, y3_left, y3_right, tol)
    return quadgk(
        y2 -> mesure_integral_y3(y1, y2, x1, x2, x3, y3_left, y3_right, tol)[1],
        y2_left,
        y2_right;
        rtol=tol,
    )
end

function _mesure_frac_quadgk(x1, x2, x3, tol, bound)
    y1_left, y1_right, y2_left, y2_right, y3_left, y3_right = mesure_bounds(bound)
    result, error = quadgk(
        y1 -> mesure_integral_y2(y1, x1, x2, x3, y2_left, y2_right, y3_left, y3_right, tol)[1],
        y1_left,
        y1_right;
        rtol=tol,
    )
    return result, error
end

function mesure_frac_quadgk(x1, x2, x3; tol=tol, bound=BIG_CONST)
    return _mesure_frac_quadgk(x1, x2, x3, tol, bound)
end

function _mesure_frac_hcubature(x1, x2, x3, tol, bound)
    y1_left, y1_right, y2_left, y2_right, y3_left, y3_right = mesure_bounds(bound)
    lower = [y1_left, y2_left, y3_left]
    upper = [y1_right, y2_right, y3_right]

    integrand(y) = mesure(y[1], y[2], y[3], x1, x2, x3)
    result, error = hcubature(integrand, lower, upper; rtol=tol)
    return result, error
end

function mesure_frac_hcubature(x1, x2, x3; tol=tol, bound=BIG_CONST)
    return _mesure_frac_hcubature(x1, x2, x3, tol, bound)
end

function mesure_frac(x1, x2, x3; tol=tol, bound=BIG_CONST, method=:hcubature)
    if method == :hcubature
        return _mesure_frac_hcubature(x1, x2, x3, tol, bound)
    elseif method == :quadgk
        return _mesure_frac_quadgk(x1, x2, x3, tol, bound)
    end

    throw(ArgumentError("Unknown method=$(method). Use :hcubature or :quadgk."))
end

function compare_mesure_methods(x1, x2, x3; tol=tol, bound=BIG_CONST)
    quadgk_elapsed = @elapsed quadgk_value = _mesure_frac_quadgk(x1, x2, x3, tol, bound)
    hcubature_elapsed = @elapsed hcubature_value = _mesure_frac_hcubature(x1, x2, x3, tol, bound)

    return (
        quadgk=(value=quadgk_value[1], error=quadgk_value[2], elapsed=quadgk_elapsed),
        hcubature=(value=hcubature_value[1], error=hcubature_value[2], elapsed=hcubature_elapsed),
    )
end

mesure_frac(1.06, 2.09, 3.12)

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    println("mesure_int.jl loaded.")
    println("Call mesure_frac_hcubature(x1, x2, x3; tol=..., bound=...) for the faster 3D integral.")
    println("Call mesure_frac_quadgk(x1, x2, x3; tol=..., bound=...) for the nested QuadGK version.")
    println("Call compare_mesure_methods(x1, x2, x3; tol=..., bound=...) for a direct comparison.")
end
