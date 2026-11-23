# All plotted values are evaluated by Julia. The bundled Plotly renderer
# draws the same data on screen and when exporting a vector figure.
const _plot_direct = "#0a5b82"
const _plot_reciprocal = "#ae481f"

# The specification consists of finite numeric arrays and trusted UI strings.
function _plot_json(x)
    if x === nothing
        return "null"
    elseif x isa AbstractString
        return "\"" * replace(x, "\\"=>"\\\\", "\""=>"\\\"", "\n"=>"\\n",
                              "\r"=>"\\r", "\t"=>"\\t", "<"=>"\\u003c") * "\""
    elseif x isa AbstractDict
        return "{" * join((_plot_json(k)*":"*_plot_json(v) for (k,v) in x), ",") * "}"
    elseif x isa AbstractArray || x isa Tuple
        return "[" * join(_plot_json.(x), ",") * "]"
    elseif x isa Bool || (x isa Real && isfinite(x))
        return string(x)
    end
    throw(ArgumentError("unsupported value in plot specification"))
end

function _trace(x, y, reciprocal, axis, showlegend, sparse)
    name = reciprocal ? "Обратная r*" : "Основная r"
    suffix = axis==1 ? "" : string(axis)
    Dict("type"=>"scatter", "mode"=>sparse ? "markers" : "lines",
         "marker"=>Dict("size"=>2.5), "x"=>x, "y"=>y,
         "xaxis"=>"x"*suffix, "yaxis"=>"y"*suffix, "name"=>name,
         "legendgroup"=>name, "showlegend"=>showlegend,
         "line"=>Dict("color"=>reciprocal ? _plot_reciprocal : _plot_direct,
                      "width"=>2.4, "dash"=>reciprocal ? "dot" : "solid"),
         "hovertemplate"=>"x = %{x:.8g}<br>%{y:.8g}<extra>"*name*"</extra>")
end

source_genus(::ZolotarevFourBands) = 1
source_genus(::PaperFilter) = 2
source_genus(::FilterModel) = 3
source_genus(::FoldedFourBands) = 4
source_genus(::Union{PickFourBands,CubicFourBands}) = 5

_plot_amplitude(model)=elliptic_constants(model.T)
_plot_counts(f::FilterModel) = (f.n,f.m1,f.m2-f.m1,f.n-f.m2)
_plot_target(model,b) = b==1 ? -1 : 1
_plot_description(f::FilterModel) = "Род 3 · уровни разрезов $(f.m1) и $(f.m2)"

function _figure(model; samples=241)
    samples >= 16 || throw(ArgumentError("at least 16 samples per band are required"))
    e = model.endpoints
    bands = length(e)÷2
    traces, shapes, annotations = Dict[], Dict[], Dict[]
    amplitude = _plot_amplitude(model)
    limit = 1.15 + 1.5*amplitude.mu
    μtext = amplitude.mu>0 ? @sprintf("μ ≈ %.5g",amplitude.mu) : @sprintf("log μ ≈ %.5g",amplitude.log_mu)
    layout = Dict{String,Any}("paper_bgcolor"=>"#fcfbf8","plot_bgcolor"=>"#ffffff",
        "font"=>Dict("family"=>"DejaVu Sans, Arial, sans-serif","size"=>13,"color"=>"#253440"),
        "height"=>680,"margin"=>Dict("l"=>76,"r"=>32,"t"=>76,"b"=>62),
        "legend"=>Dict("orientation"=>"h","x"=>1,"xanchor"=>"right","y"=>1.15),
        "title"=>Dict("text"=>"Род $(source_genus(model)) · n = $(model.n) · $μtext","x"=>0.06,"y"=>0.98),
        "hovermode"=>"closest","shapes"=>shapes,"annotations"=>annotations)
    axis = Dict("showline"=>true,"linecolor"=>"#9aa5aa","gridcolor"=>"#e6e9e9",
                "zeroline"=>false,"ticks"=>"outside","ticklen"=>3,"automargin"=>true)
    layout["xaxis"]=merge(axis,Dict("domain"=>[0,1],"anchor"=>"y","title"=>"x","nticks"=>8))
    layout["yaxis"]=merge(axis,Dict("domain"=>[.54,1],"anchor"=>"x","title"=>"Значение дроби","range"=>[-limit,limit],"nticks"=>7))
    counts = _plot_counts(model)
    for b in 1:bands
        a,z = e[2b-1],e[2b]
        # Cosine spacing resolves the square-root coordinate near band ends.
        x=[clamp(a+(z-a)*sinpi(j/(2*(samples-1)))^2,a,z) for j in 0:samples-1]
        x[1]=a; x[end]=z
        rows = [wall_values(model, t,b) for t in x]
        sparse = counts[b] > (samples-1)÷8
        target = _plot_target(model,b)
        push!(shapes,Dict("type"=>"rect","xref"=>"x","yref"=>"paper",
              "x0"=>a,"x1"=>z,"y0"=>.54,"y1"=>1,"fillcolor"=>b==1 ? "#eef1f5" : "#e3efe9",
              "layer"=>"below","line"=>Dict("width"=>0)))
        push!(shapes,Dict("type"=>"line","x0"=>a,"x1"=>z,"y0"=>target,"y1"=>target,
                   "layer"=>"below","line"=>Dict("color"=>"#818885","dash"=>"dash","width"=>1)))
        width = (1-.06*(bands-1))/bands
        left = (b-1)*(width+.06)
        title = b==1 ? "Подавление" : "Пропускание $(b-1)"
        push!(annotations, Dict("text"=>title,"xref"=>"paper","yref"=>"paper",
              "x"=>left+width/2,"y"=>.35,"xanchor"=>"center","showarrow"=>false))
        k = string(b+1)
        layout["xaxis"*k]=merge(axis,Dict("domain"=>[left,left+width],"anchor"=>"y"*k,
            "tickvals"=>[a,z],"ticktext"=>[@sprintf("%.5g",a),@sprintf("%.5g",z)]))
        layout["yaxis"*k]=merge(axis,Dict("domain"=>[.02,.29],"anchor"=>"x"*k,
            "range"=>[-1.12,1.12],"tickvals"=>[-1,0,1],"showticklabels"=>b==1))
        for reciprocal in (false,true)
            v = reciprocal ? :reciprocal : :direct
            err = reciprocal ? :reciprocal_ripple : :ripple
            push!(traces,_trace(x,real.(getproperty.(rows,v)),reciprocal,1,b==1,sparse))
            push!(traces,_trace(x,getproperty.(rows,err),reciprocal,b+1,false,sparse))
        end
    end
    # Omit values beyond the viewport at real poles so the plot does not
    # connect the two sides of an asymptote.
    for gap in 1:bands-1
        a,z=e[2gap],e[2gap+1]
        xs=[clamp(a+(z-a)*sinpi(j/(2*(samples-1)))^2,a,z) for j in 0:samples-1]
        rows=map(xs) do x
            try
                evaluate(model,x)
            catch err
                err isa DomainError || rethrow()
                nothing
            end
        end
        for reciprocal in (false,true)
            key=reciprocal ? :reciprocal : :direct
            ys=map(rows) do row
                row===nothing && return nothing
                v=real(getproperty(row,key))
                abs(v)<=1.4*limit ? v : nothing
            end
            trace=_trace(xs,ys,reciprocal,1,false,false)
            trace["line"]["width"]=1.4
            trace["opacity"]=0.78
            push!(traces,trace)
        end
    end
    push!(annotations,Dict("text"=>"Ошибка относительно цели, в единицах μ",
        "xref"=>"paper","yref"=>"paper","x"=>0,"y"=>.411,"xanchor"=>"left","showarrow"=>false))
    Dict("data"=>traces,"layout"=>layout)
end

"""Write an offline, interactive HTML plot with SVG export and four ripple panels."""
function plot_filter(path, model; samples=241)
    spec=_plot_json(_figure(model;samples))
    renderer=read(joinpath(@__DIR__, "../../vendor/plotly/plotly-basic.min.js"),String)
    c=_plot_amplitude(model)
    amplitude = c.mu>0 ? @sprintf("μ = %.6g",c.mu) : @sprintf("log μ = %.6g",c.log_mu)
    open(path,"w") do io
        print(io,"""<!doctype html><html lang="ru"><meta charset="utf-8">
        <meta name="viewport" content="width=device-width,initial-scale=1">
        <title>Фильтр · степень $(model.n)</title>
        <style>body{margin:0;background:#fcfbf8;color:#253440;font:16px system-ui,sans-serif}
        main{max-width:1200px;margin:25px auto;padding:0 16px}h1{font-size:clamp(22px,3vw,34px);margin:0}
        p{margin:8px 0 16px}#chart{width:100%;min-width:680px}.panels{overflow-x:auto}</style>
        <main><h1>Фильтр степени $(model.n) · полосы: $(length(model.endpoints)÷2)</h1>
        <p> $(_plot_description(model)) · Im τ = $(@sprintf("%.7g",model.T)) · $amplitude</p>
        <div class="panels"><div id="chart" role="img" aria-label="Две дроби и их нормированные ошибки на рабочих полосах"></div></div>
        <p>Цветные области сверху обозначают рабочие полосы в общем масштабе x. Внизу каждый интервал увеличен отдельно.
        Для полос с плотными колебаниями показаны отсчёты. Выделите область для увеличения.</p></main>
        <script>""",renderer,"""</script><script id="figure" type="application/json">""",spec,
        """</script><script>const f=JSON.parse(document.getElementById('figure').textContent);
        Plotly.newPlot('chart',f.data,f.layout,{responsive:true,displaylogo:false,
        toImageButtonOptions:{format:'svg',filename:'four-band-filter'}});</script></html>""")
    end
    path
end
