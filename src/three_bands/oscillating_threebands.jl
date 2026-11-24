"""
    oscillating_threebands(n, epsilon; order=28)

Construct three free real bands with targets (-1, 1, 1) and weights (1, 1, 1/2).
Compare the two seam phases and place the two inner cuts using their common
contact window. The interface and the endpoints are computed in Float64.
Use `welded_value(solution, x)` on any of the three closed bands.
"""
function oscillating_threebands(n::Integer, epsilon::Real; order::Int=28)
    m = 2 * fld(n, 4)
    n >= 4 || throw(ArgumentError("at least four tiles are needed for the two blocks"))
    phases = map((0, 1)) do offset
        welded_threebands(n, m, epsilon; third_weight=0.5, offset, order)
    end
    chosen = phases[1].modulus >= phases[2].modulus ? phases[1] : phases[2]
    outer = elliptic_constants(2 * chosen.modulus)
    outer.k >= sqrt(floatmin(Float64)) ||
        throw(ArgumentError("the symmetric outer coordinates require higher precision"))

    marks = map(phases) do f
        previous = welded_wall(f, 1, 1, m - 1).x
        narrow = welded_wall(f, 1, 1, m).x
        broad = welded_wall(f, 2, 1, m).x
        (; previous, narrow, broad)
    end
    L = minimum(mark.narrow for mark in marks)
    tip = minimum(mark.broad for mark in marks)
    R = L + (tip - L) / 2
    maximum(mark.previous for mark in marks) < L < R < tip ||
        throw(ArgumentError("the phases have no resolved common contact window"))

    grids = (left=_weld_grid(chosen, :left), right=_weld_grid(chosen, :right),
             tip=_weld_grid(chosen, :tip))
    filter = merge(chosen, (; grids))
    target = (M=chosen.modulus, q=outer.k, K=outer.K)
    ends = [-1.0, -outer.k, outer.k, L, R, 1.0]
    (; filter, target, ends, marks, moduli=map(f -> f.modulus, phases))
end
