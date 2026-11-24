include("../src/four_bands/StiefelFourBands.jl")
using .StiefelFourBands
using Printf

n = isempty(ARGS) ? 16 : parse(Int, ARGS[1])
epsilon = length(ARGS) < 2 ? 0.1 : parse(Float64, ARGS[2])
weights = length(ARGS) <= 3 ? (1.0, 1.0, 0.5) : parse.(Float64, ARGS[4:end])
solution = oscillating_threebands(n, epsilon; weights)
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

include("three_bands/samples.jl")

traces = []
values = []
errors = []
view = isnothing(periods) ? "полные полосы" : "до $periods периодов на полосу"
layout = Dict{String,Any}("height" => 750, "showlegend" => false,
    "title" => "n = $n · ε = $epsilon · $view",
    "grid" => Dict("rows" => 3, "columns" => 1, "pattern" => "independent"))
for band in 1:3
    samples = sample_threeband(solution, band, periods)
    logx = samples.logx
    push!(values, samples.value)
    push!(errors, samples.error)
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
        width = epsilon / solution.weights[band]
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
