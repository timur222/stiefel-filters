#!/usr/bin/env julia
include("../src/four_bands/StiefelFourBands.jl")
using .StiefelFourBands
using Printf

const USAGE = """
Julia four-band filter
  julia --project=. bin/fourbands.jl demo [options]
  julia --project=. bin/fourbands.jl solve N E0 E1 E2 E3 E4 E5 E6 E7 [options]
  julia --project=. bin/fourbands.jl load filter.toml [options]

  julia --project=. bin/fourbands.jl periods periods.toml [options]

  --plot PATH         Offline interactive HTML (default four-bands.html)
  --save PATH         Save compact filter parameters as TOML
  --at RE IM          Evaluate both fractions, repeatable
  --method METHOD     sc or theta for genus 3, sc for genus 2, elliptic for genus 1, rational for genera 4 and 5
  --genus 1|2|3|5        Family for solve with prescribed endpoints (default 3)
  --levels M1 M2      Fix the slit levels during synthesis
  --sides S1 S2       Fix trim sides, above or below
  --order N           Quadrature order for synthesis (default 128)
  --samples N         Number of visual samples per band (default 241)
"""

function run(args)
    if isempty(args) || first(args) in ("-h","--help")
        print(USAGE)
        return
    end
    input=copy(args)
    poparg()=isempty(input) ? throw(ArgumentError("missing argument; use --help")) : popfirst!(input)
    command=poparg()
    command in ("demo","solve","load","periods") || throw(ArgumentError("expected demo, solve, load, or periods"))
    degree=command=="solve" ? parse(Int,poparg()) : 3
    endpoints=command=="solve" ? [parse(Float64,poparg()) for _ in 1:8] :
                                  [-1.,0.,.10525,.10772,.35397,.37827,.94801,1.]
    source=command in ("load","periods") ? poparg() : nothing
    plotpath="four-bands.html"; savepath=nothing; levels=nothing; sides=nothing
    genus=3; order=128; method=nothing; points=ComplexF64[]; samples=241
    while !isempty(input)
        opt=poparg()
        command=="periods" && opt in ("--sides","--levels") && throw(ArgumentError("set slit levels and sides in the period file"))
        if command=="load" && opt in ("--levels","--sides","--order")
            throw(ArgumentError("$opt applies to synthesis only"))
        end
        if opt=="--plot"
            plotpath=poparg()
        elseif opt=="--save"
            savepath=poparg()
        elseif opt=="--at"
            push!(points,complex(parse(Float64,poparg()),parse(Float64,poparg())))
        elseif opt=="--method"
            method=Symbol(poparg())
            method in (:theta,:sc,:elliptic,:rational) || throw(ArgumentError("expected --method theta, sc, elliptic or rational"))
        elseif opt=="--sides"
            sides=(Symbol(poparg()),Symbol(poparg()))
        elseif opt=="--genus"
            command=="solve" || throw(ArgumentError("--genus applies to solve; period files specify their own genus"))
            genus=parse(Int,poparg())
            genus in (1,2,3,5) || throw(ArgumentError("supported genera: 1, 2, 3 and 5"))
        elseif opt=="--levels"
            levels=(parse(Int,poparg()),parse(Int,poparg()))
        elseif opt=="--order"
            order=parse(Int,poparg())
        elseif opt=="--samples"
            samples=parse(Int,poparg())
        else
            throw(ArgumentError("unknown option $opt"))
        end
    end
    samples >= 16 || throw(ArgumentError("--samples must be at least 16"))
    if savepath !== nothing && abspath(savepath)==abspath(plotpath)
        throw(ArgumentError("use different paths for the filter and the plot"))
    end
    source !== nothing && abspath(source)==abspath(plotpath) &&
        throw(ArgumentError("plot path overlaps saved filter"))
    genus!=3 && (levels!==nothing || sides!==nothing) && throw(ArgumentError("--levels and --sides apply to the genus-three inverse"))
    model = if command=="periods"
        read_period_filter(source;order)
    elseif command=="load"
        load_filter(source)
    elseif genus==5
        degree==3 || throw(ArgumentError("the explicit genus-five inverse requires degree 3"))
        cubic_fourbands(endpoints)
    elseif genus==1
        zolotarev_fourbands(endpoints,degree)
    elseif genus==2
        paper_fourbands(endpoints,degree;order)
    else
        synthesize(endpoints,degree;levels,sides,order)
    end
    c=StiefelFourBands._plot_amplitude(model)
    if model isa FilterModel
        println("genus=3   m=(",model.m1,", ",model.m2,")   trims=",model.trims,"   sides=",model.sides)
        command=="periods" && println("Ω = ",period_parameters(model).omega)
        println("checks: ",model.diagnostics)
    elseif model isa PaperFilter
        println("genus=2   m=",model.m,"   Ω = ",model.omega)
    elseif model isa ZolotarevFourBands
        println("genus=1   Ω₁₁ = ",model.H/2)
    elseif model isa FoldedFourBands
        println("genus=4   cubic rational function")
    else
        println("genus=5   cubic rational function")
    end
    @printf("n=%d   τ=%.13gi   μ=%.13g   log μ=%.13g\n",model.n,model.T,c.mu,c.log_mu)
    println("endpoints: ",model.endpoints)
    if isnothing(method)
        method = model isa ZolotarevFourBands ? :elliptic : model isa Union{PaperFilter,FilterModel} ? :sc : :rational
    end
    methods=model isa FilterModel ? (:theta,:sc) : model isa PaperFilter ? (:sc,) : model isa ZolotarevFourBands ? (:elliptic,) : (:rational,)
    method in methods || throw(ArgumentError("evaluation methods for this genus: $methods"))
    savepath!==nothing && !(model isa FilterModel) && throw(ArgumentError("compact --save is for genus 3; reuse the input TOML or endpoint arguments for this genus"))
    source!==nothing && savepath!==nothing && command=="periods" && abspath(source)==abspath(savepath) &&
        throw(ArgumentError("saved filter would overwrite the input periods"))
    curve = method==:theta && any(z->!isreal(z),points) ? ThetaCurve(model.geom) : nothing
    for z in points
        v=method==:theta ? evaluate(model,z;method,curve) : evaluate(model,z;method)
        println("x=",z,"   r=",v.direct,"   r*=",v.reciprocal)
    end
    savepath===nothing || println("TOML: ",save_filter(savepath,model))
    println("HTML: ",plot_filter(plotpath,model;samples))
end

try
    run(ARGS)
catch ex
    ex isa InterruptException && rethrow()
    println(stderr,sprint(showerror,ex))
    exit(1)
end
