"""Single-slit genus-two filter, optionally restricted to four working intervals."""
struct PaperFilter
    n::Int
    m::Int
    T::Float64
    omega::Matrix{Float64}
    geom::Geometry
    endpoints::Vector{Float64}
    walls::Vector{Int}
end

function _paper_geometry(e,order)
    M=reduce(vcat, [transpose(moments(e,i,order)) for i in 1:5])
    p=[M[[2,4],1:2] \ [2.,0.];0.]
    H=[dot(M[3,:],p),-dot(M[5,:],p)]
    left=dot(M[1,:],p)
    root=-p[1]/p[2]
    e[4]<root<e[5] || throw(GeometryError("genus-two slit tip outside its gap"))
    slit=1-dot(p,moments(e,4,order,_angle(e[4],e[5],root),π))
    -1<slit<1 || throw(GeometryError("degenerate genus-two slit"))
    Geometry(e,p,H,[root],[slit],left,order,M)
end

function _paper_height(g, m, T, x, wall)
    i = 2wall-1
    a,b = g.e[i:i+1]
    a <= x <= b || throw(ArgumentError("point outside the selected genus-two wall"))
    increment = dot(g.p,moments(g.e,i,g.order,_angle(a,b,x),π))
    wall==1 && return g.left-increment
    wall==2 && return increment
    m*T-increment
end

function _paper_cut(g,m,T,wall,tiles)
    a,b = g.e[2wall-1:2wall]
    for _ in 1:56
        x = (a+b)/2
        if _paper_height(g,m,T,x,wall) < tiles*T
            a=x
        else
            b=x
        end
    end
    (a+b)/2
end

"""
    paper_filter(tau, n, m, omega22; split=nothing, order=128)

Genus-two single-slit solution from (17), (23). We plot the finite projective
coordinate (x−1)/(x+1) of the paper. To retain four bands use
`split=(passband, tile, from, to)`: remove the open interval between heights
(tile+from)T and (tile+to)T on passband 1 or 2, strictly inside one tile.
"""
function paper_filter(tau,n::Integer,m::Integer,omega22::Real;split=nothing,order::Int=128)
    Pi=paper_period_matrix(tau,n,m,omega22)
    T=2imag(Pi[1,1])/n
    X=paper_branchpoints(Pi)
    e=[[(x-1)/(x+1) for x in X[1:5]];1.]
    all(>(0),diff(e)) || throw(GeometryError("genus-two branch points coalesce in this chart"))
    g=_paper_geometry(e,order)
    measured=[g.left;g.H]/T
    maximum(abs.(measured-[n,m,n-m]))<1e-7 ||
        throw(GeometryError("genus-two height check failed in tile units; inspect conditioning"))
    E=copy(e)
    walls=[1,2,3]
    if !isnothing(split)
        length(split)==4 || throw(ArgumentError("split requires passband, tile, from, to"))
        pass, tile, lo, hi = split
        pass isa Integer && pass in (1,2) && tile isa Integer && 0<lo<hi<1 ||
            throw(ArgumentError("choose passband 1 or 2, integer tile, and 0 < from < to < 1"))
        start,stop = pass==1 ? (0,m) : (m,n)
        start<=tile<stop || throw(ArgumentError("tile outside the requested passband"))
        wall=pass+1
        a=_paper_cut(g,m,T,wall,tile+lo)
        b=_paper_cut(g,m,T,wall,tile+hi)
        e[2wall-1]<a<b<e[2wall] || throw(GeometryError("split endpoints outside available precision"))
        E=[e[1:2wall-1];a;b;e[2wall:end]]
        insert!(walls,wall,wall)
    end
    PaperFilter(n,m,T,imag.(Pi),g,E,walls)
end

function height_tiles(f::PaperFilter,x::Real,band::Integer)
    _paper_height(f.geom,f.m,f.T,x,f.walls[band])/f.T
end

function evaluate(f::PaperFilter,z::Number;method=:sc,curve=nothing)
    method==:sc || throw(ArgumentError("genus-two periods and endpoints use theta; point evaluation uses :sc"))
    isfinite(z) || throw(ArgumentError("evaluation requires a finite point"))
    g=f.geom
    if isreal(z)
        x=real(z)
        for w in 1:3
            if g.e[2w-1]<=x<=g.e[2w]
                t=_paper_height(g,f.m,f.T,x,w)/f.T
                return _wall_values_at(t,f.T,w==1 ? -1 : 1)
            end
        end
        u = if g.e[2]<x<g.e[3] || g.e[4]<x<g.e[5]
            onbase = x<g.e[3]
            i=onbase ? 2 : 4
            inc=dot(g.p,moments(g.e,i,g.order,_angle(g.e[i],g.e[i+1],x),π))
            complex(onbase ? -1+inc : 1-inc,onbase ? 0. : f.m*f.T)
        else
            complex(_outer_abscissa(g,x),f.n*f.T)
        end
        level = x < g.e[1] || x > g.e[end] ? f.n : x < g.e[3] ? 0 : f.m
        abs(real(u)) > 2e-12 || throw(DomainError(x,"pole or unresolved real SC coordinate"))
        a,b=_pair(real(u),f.T,elliptic_constants(f.T))
        direct, reciprocal = isodd(level) ? (real(b),real(a)) : (real(a),real(b))
        return (;direct,reciprocal,method)
    else
        u=_upper_coordinate(g,complex(real(z),abs(imag(z))))
    end
    direct,reciprocal=_pair(u,f.T,elliptic_constants(f.T))
    imag(z)<0 && return (;direct=conj(direct),reciprocal=conj(reciprocal),method)
    (;direct,reciprocal,method)
end

_plot_counts(f::PaperFilter) = [w==1 ? f.n : w==2 ? f.m : f.n-f.m for w in f.walls]
_plot_description(f::PaperFilter) = "Род 2 · уровень разреза $(f.m)"

"""
    paper_filter(tau, P; kwargs...)

Read the complete marked genus-two Ω or Π=iΩ. Recover n,m from (17),
then use the same projection and period checks as the four-argument input.
"""
function paper_filter(tau,P::AbstractMatrix;kwargs...)
    size(P)==(2,2) && all(x->x isa Number && isfinite(x),P) ||
        throw(ArgumentError("expected a finite 2×2 marked matrix"))
    omega = if all(isreal,P)
        real.(P)
    elseif all(z->iszero(real(z)),P)
        imag.(P)
    else
        throw(ArgumentError("expected real Ω or purely imaginary Π"))
    end
    maximum(abs,omega-transpose(omega))<=5e-13*max(1.,maximum(abs,omega)) ||
        throw(ArgumentError("Ω must be symmetric"))
    T=_period_tau(tau)
    tiles=2 .* omega[1,:] ./ T
    all(x->0<x<2^53,tiles) || throw(ArgumentError("invalid heights in first row"))
    n,k=round.(Int,tiles)
    maximum(abs,tiles-[n,k])<=1e-7 || throw(ArgumentError("genus-two periods do not satisfy (17) with integer n,m"))
    paper_filter(im*T,n,n-k,omega[2,2];kwargs...)
end
