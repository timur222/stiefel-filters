"""Write a reusable compact filter in human-readable TOML."""
function save_filter(path, model::FilterModel)
    p = elliptic_constants(model.T)
    data = Dict("format_version"=>1, "n"=>model.n, "m1"=>model.m1, "m2"=>model.m2,
                "endpoints"=>model.endpoints, "e"=>model.geom.e, "order"=>model.geom.order,
                "T"=>model.T, "sides"=>string.(collect(model.sides)), "trims"=>collect(model.trims),
                "p"=>model.geom.p, "slit"=>model.geom.slit, "mu"=>p.mu, "log_mu"=>p.log_mu)
    open(io->TOML.print(io, data; sorted=true), path, "w")
    path
end

"""Rebuild a saved Julia filter and recheck its periods and trims."""
function load_filter(path)
    d = TOML.parsefile(path)
    get(d,"format_version",0)==1 || throw(ArgumentError("unrecognized filter format"))
    restore_model(; e=d["e"], endpoints=d["endpoints"], n=d["n"], m1=d["m1"], m2=d["m2"],
                  sides=d["sides"], trims=d["trims"], T=d["T"], order=d["order"])
end
