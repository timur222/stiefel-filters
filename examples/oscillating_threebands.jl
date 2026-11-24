include("../src/four_bands/StiefelFourBands.jl")
using .StiefelFourBands

n = isempty(ARGS) ? 16 : parse(Int, ARGS[1])
epsilon = length(ARGS) < 2 ? 0.1 : parse(Float64, ARGS[2])
solution = oscillating_threebands(n, epsilon)
println("Концы полос = ", solution.ends)
println("Модули двух склеек = ", solution.moduli)

# Separate axes show the ripples even when the physical endpoints are tiny.
traces = []
layout = Dict{String,Any}("height" => 750, "showlegend" => false,
    "title" => "Взвешенная ошибка / ε · n = $n · ε = $epsilon",
    "grid" => Dict("rows" => 3, "columns" => 1, "pattern" => "independent"))
for band in 1:3
    a, b = solution.ends[2band-1:2band]
    logx = range(log10(abs(a)), log10(abs(b)); length=max(300, 20n))
    xs = clamp.(sign(a) .* (10.0 .^ logx), a, b)
    xs[1], xs[end] = a, b
    ys = [welded_value(solution, x).normalized for x in xs]
    suffix = band == 1 ? "" : string(band)
    push!(traces, Dict("type" => "scatter", "mode" => "lines", "name" => "E$band",
        "xaxis" => "x$suffix", "yaxis" => "y$suffix", "x" => collect(logx), "y" => ys))
    layout["xaxis$suffix"] = Dict("title" => "log₁₀|x|", "range" => [first(logx), last(logx)])
    layout["yaxis$suffix"] = Dict("title" => "E$band", "range" => [-1.1, 1.1],
                                  "tickvals" => [-1, 0, 1])
end

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
