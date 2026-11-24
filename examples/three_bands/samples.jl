# A full period spans two successive tiles. Select physical endpoints of
# the window using tile numbers, without visiting the preceding extrema.
function plot_nodes(solution, band, periods)
    f = solution.filter
    chart_band = solution.reflected && band != 1 ? 5-band : band
    lo, hi = ((0, f.n), (0, f.m), (f.m, f.n))[chart_band]
    function abscissa(j)
        block = chart_band == 3 || j > f.m ? 2 : 1
        wall = band == 1 ? -1 : 1
        x = StiefelFourBands.welded_wall(f, block, wall, j, solution.target).x
        solution.reflected ? solution.target.q / x : x
    end
    if !isnothing(periods) && periods < cld(hi - lo, 2)
        lo += fld(hi - lo - 2periods, 2)
        return sort!(abscissa.(lo:lo+2periods))
    end
    a, b = solution.ends[2band-1:2band]
    contacts = filter(x -> a < x < b, abscissa.(lo:hi))
    sort!([a; contacts; b])
end

function sample_threeband(solution, band, periods=nothing)
    nodes = plot_nodes(solution, band, periods)
    steps = max(20, cld(300, length(nodes)-1))
    xs = [first(nodes)]
    for (a,b) in zip(nodes, Iterators.drop(nodes, 1))
        # Keep the exact physical contacts in the data shown by the browser.
        u, v = log(abs(a)), log(abs(b))
        for k in 1:steps-1
            push!(xs, clamp(sign(a)*exp(u+(v-u)*(k/steps)), a, b))
        end
        push!(xs, b)
    end
    logx = Float64.(log10.(abs.(xs)))
    rows = [welded_value(solution, x) for x in xs]
    (; x=xs, logx, value=[r.value for r in rows], error=[r.normalized for r in rows])
end
