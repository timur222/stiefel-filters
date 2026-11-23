function _elliptic_period_degree(tau,P)
    size(P)==(1,1) || throw(ArgumentError("expected a 1×1 scalar period matrix"))
    z=only(P)
    z isa Number && isfinite(z) && (isreal(z) || iszero(real(z))) ||
        throw(ArgumentError("expected real Ω or purely imaginary Π"))
    h=2*(isreal(z) ? real(z) : imag(z))/_period_tau(tau)
    0<h<2^53 && abs(h-round(h))<=1e-7 || throw(ArgumentError("2Ω₁₁/Im τ must be a positive integer"))
    round(Int,h)
end

"""
    filter_from_periods(tau, P; kwargs...)

Construct a filter in the implemented marked family of genus size(P,1)
for real Ω or purely imaginary Π=iΩ, given the independent elliptic τ.
Use gaps for genus 1, split for genus 2, or trims/sides for genus 3.
"""
function filter_from_periods(tau,P::AbstractMatrix;kwargs...)
    g=size(P,1)
    size(P,2)==g || throw(ArgumentError("period matrix must be square"))
    g==1 && return zolotarev_fourbands(im*_period_tau(tau),_elliptic_period_degree(tau,P);kwargs...)
    g==2 && return paper_filter(tau,P;kwargs...)
    g==3 && return _three_from_periods(tau,P;kwargs...)
    g==4 && return folded_from_periods(tau,P;kwargs...)
    g==5 && return pick_from_periods(tau,P;kwargs...)
    throw(ArgumentError("supported marked families have genus 1 to 5"))
end

"""Read a TOML specification of the elliptic period and the marked source curve."""
function read_period_filter(path; order=128)
    d=TOML.parsefile(path)
    g=get(d,"genus",nothing)
    g isa Integer && !(g isa Bool) && g in 1:5 || throw(ArgumentError("genus must be between 1 and 5"))
    keys_for_mode = if g==1
        haskey(d,"omega") ? ["omega","gaps"] : ["n","gaps"]
    elseif g==2
        haskey(d,"omega") ? ["omega","split"] : ["n","m","omega22","split"]
    elseif g==4
        haskey(d,"omega") ? ["omega"] : ["scale","turning","odd_root"]
    elseif g==5
        haskey(d,"omega") ? ["omega"] : ["poles","residues","offset"]
    elseif haskey(d,"omega")
        ["omega","trims","sides"]
    else
        ["n","m1","m2","omega22","omega23","omega33","trims","sides"]
    end
    extra=setdiff(keys(d),["genus";"tau_im";keys_for_mode])
    isempty(extra) || throw(ArgumentError("unexpected fields for genus $g: $(join(extra,", "))"))
    T=get(d,"tau_im",nothing)
    T isa Real && !(T isa Bool) && isfinite(T) && T>0 ||
        throw(ArgumentError("tau_im must be positive (τ = i*tau_im)"))
    tau=im*T
    require(key)=haskey(d,key) ? d[key] : throw(ArgumentError("missing field $key"))
    function matrix()
        rows=require("omega")
        rows isa AbstractVector && length(rows)==g &&
            all(row->row isa Vector && length(row)==g && all(v->v isa Real && isfinite(v),row),rows) ||
            throw(ArgumentError("omega must have $g rows of $g finite real entries"))
        reduce(vcat,permutedims.(rows))
    end
    if g==1
        n=haskey(d,"omega") ? _elliptic_period_degree(tau,matrix()) : require("n")
        return haskey(d,"gaps") ? zolotarev_fourbands(tau,n;gaps=d["gaps"]) : zolotarev_fourbands(tau,n)
    elseif g==4
        return haskey(d,"omega") ? folded_from_periods(tau,matrix();order) :
            folded_fourbands(tau;scale=get(d,"scale",1.),turning=get(d,"turning",2.),odd_root=get(d,"odd_root",-3.))
    elseif g==5
        if haskey(d,"omega")
            return pick_from_periods(tau,matrix();order)
        end
        return pick_fourbands(tau;poles=get(d,"poles",(-1.,1.)),
            residues=get(d,"residues",(1.,1.)),offset=get(d,"offset",0.))
    elseif g==2
        split = get(d,"split",nothing)
        if !isnothing(split)
            cols = ["passband","tile","from","to"]
            split isa AbstractDict && Set(keys(split))==Set(cols) ||
                throw(ArgumentError("split requires passband, tile, from and to"))
            split=Tuple(split[k] for k in cols)
        end
        haskey(d,"omega") && return paper_filter(tau,matrix();split,order)
        return paper_filter(tau,require("n"),require("m"),require("omega22");split,order)
    end
    trims=get(d,"trims",(0.,0.))
    sides=get(d,"sides",(:above,:above))
    if haskey(d,"omega")
        return filter_from_periods(tau,matrix();trims,sides,order)
    end
    filter_from_periods(tau,require("n"),require("m1"),require("m2");
        omega22=require("omega22"),omega33=require("omega33"),omega23=get(d,"omega23",nothing),
        trims,sides,order)
end
