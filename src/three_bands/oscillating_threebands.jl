"""
    oscillating_threebands(n, epsilon; order=28)

Construct three free real bands with targets (-1, 1, 1) and weights (1, 1, 1/2).
Compare the two seam phases and place the two inner cuts using their common
contact window. The interface uses Float64. Long rectangles use BigFloat
for the physical endpoints so that distinct bands remain representable.
Use `welded_value(solution, x)` on any of the three closed bands.
"""
function oscillating_threebands(n::Integer, epsilon::Real; order::Int=28)
    m = 2 * fld(n, 4)
    n >= 4 || throw(ArgumentError("at least four tiles are needed for the two blocks"))
    phases = map((0, 1)) do offset
        welded_threebands(n, m, epsilon; third_weight=0.5, offset, order)
    end
    # Both phases share the same two bulk heights, including for odd n.
    selected = phases[1].excess >= phases[2].excess ? 1 : 2
    chosen = phases[selected]
    targets = map(f -> _weld_outer(f.modulus), phases)
    target = targets[selected]
    marks = map(phases, targets) do f, outer
        previous = welded_wall(f, 1, 1, m - 1, outer).x
        narrow = welded_wall(f, 1, 1, m, outer).x
        broad = welded_wall(f, 2, 1, m, outer).x
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
    seam_heights(wall) = map(i -> welded_wall_height(filter, i, wall, m), (1, 2))
    seam = (left=seam_heights(-1), right=seam_heights(1))
    ends = [-1.0, -target.q, target.q, L, R, 1.0]
    (; filter, target, ends, marks, seam, moduli=map(f -> f.modulus, phases))
end
