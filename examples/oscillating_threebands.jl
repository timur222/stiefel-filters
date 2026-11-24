include("../src/four_bands/StiefelFourBands.jl")
using .StiefelFourBands
using Printf

n = isempty(ARGS) ? 16 : parse(Int, ARGS[1])
epsilon = length(ARGS) < 2 ? 0.1 : parse(Float64, ARGS[2])
solution = oscillating_threebands(n, epsilon)
periods = if length(ARGS) < 3
    solution.target.q isa BigFloat ? 6 : nothing
elseif ARGS[3] == "all"
    nothing
else
    parse(Int, ARGS[3])
end
isnothing(periods) || periods > 0 || error("Число периодов должно быть положительным")
println("Концы полос = [", join((@sprintf("%.12g", x) for x in solution.ends), ", "), "]")
println("Модули двух склеек = ", solution.moduli)

# A full period spans two successive tiles. Select physical endpoints of
# the window using tile numbers, without visiting the preceding extrema.
function plot_nodes(solution, band, periods)
    f = solution.filter
    lo, hi = ((0, f.n), (0, f.m), (f.m, f.n))[band]
    function abscissa(j)
        block = band == 3 || j > f.m ? 2 : 1
        wall = band == 1 ? -1 : 1
        StiefelFourBands.welded_wall(f, block, wall, j, solution.target).x
    end
    if !isnothing(periods) && periods < cld(hi - lo, 2)
        lo += fld(hi - lo - 2periods, 2)
        return sort!(abscissa.(lo:lo+2periods))
    end
    a, b = solution.ends[2band-1:2band]
    contacts = filter(x -> a < x < b, abscissa.(lo:hi))
    sort!([a; contacts; b])
end

traces = []
values = []
errors = []
view = isnothing(periods) ? "полные полосы" : "до $periods периодов на полосу"
layout = Dict{String,Any}("height" => 750, "showlegend" => false,
    "title" => "n = $n · ε = $epsilon · $view",
    "grid" => Dict("rows" => 3, "columns" => 1, "pattern" => "independent"))
for band in 1:3
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
    push!(values, [r.value for r in rows])
    push!(errors, [r.normalized for r in rows])
    suffix = band == 1 ? "" : string(band)
    push!(traces, Dict("type" => "scatter", "mode" => "lines", "name" => "E$band",
        "xaxis" => "x$suffix", "yaxis" => "y$suffix", "x" => logx, "y" => errors[end]))
    layout["xaxis$suffix"] = Dict("title" => "log₁₀|x|", "range" => [first(logx), last(logx)])
    layout["yaxis$suffix"] = Dict("title" => "E$band", "range" => [-1.1, 1.1],
                                  "tickvals" => [-1, 0, 1])
end

buttons = []
for (label, data, normalized) in (("Взвешенная ошибка / ε", errors, true), ("Дробь f(x)", values, false))
    axes = Dict{String,Any}()
    for band in 1:3
        axis = band == 1 ? "yaxis" : "yaxis$band"
        center = band == 1 ? -1 : 1
        width = (band == 3 ? 2 : 1) * epsilon
        axes["$axis.range"] = normalized ? [-1.1, 1.1] : center .+ [-1.1width, 1.1width]
        axes["$axis.tickvals"] = normalized ? [-1, 0, 1] : [center-width, center, center+width]
    end
    push!(buttons, Dict("label" => label, "method" => "update", "args" => (Dict("y" => data), axes)))
end
layout["updatemenus"] = [Dict("type" => "buttons", "direction" => "right",
    "x" => 1, "xanchor" => "right", "y" => 1.12, "buttons" => buttons)]

output = joinpath(mkpath(joinpath(@__DIR__, "../output")), "oscillating_threebands.html")
renderer = read(joinpath(@__DIR__, "../vendor/plotly/plotly-basic.min.js"), String)
figure = StiefelFourBands._plot_json(Dict("data" => traces, "layout" => layout))
open(output, "w") do io
    print(io, """<!doctype html><html lang="ru"><meta charset="utf-8">
        <title>Колебания на трёх полосах</title><div id="plot"></div><script>""",
        renderer, "</script><script>const fig=", figure,
        ";Plotly.newPlot('plot',fig.data,fig.layout,{responsive:true});</script></html>")
end
println("График сохранён в ", abspath(output))
