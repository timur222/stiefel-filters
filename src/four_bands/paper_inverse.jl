# In the article's homology basis a-periods are the doubled alternating
# integrals over (e₂,e₃),(e₄,e₅). The b-periods are cumulative band
# integrals. The common choice of signs gives 2Ω₁⋅ = (n,n−m) Im τ.
function _paper_omega(g)
    M = g.mom[:,1:2]
    A = 2 .* vcat(transpose(M[2,:]),-transpose(M[4,:]))
    B = 2 .* vcat(transpose(M[1,:]),transpose(M[1,:]-M[3,:]))
    raw = B / A
    maximum(abs,raw-transpose(raw)) <= 1e-9*max(1,maximum(abs,raw)) ||
        throw(GeometryError("genus-two period matrix failed its symmetry check"))
    Ω=(raw+transpose(raw))/2
    0<Ω[1,2]<min(Ω[1,1],Ω[2,2]) || throw(GeometryError("outside the marked genus-two period cone"))
    Ω
end

function _paper_gap_check(g,m,n,T,wall,a,b)
    tile = floor(Int,_paper_height(g,m,T,a,wall)/T)
    start,stop = wall==2 ? (0,m) : (m,n)
    for j in max(start,tile-1):min(stop,tile+2)
        xj = _paper_cut(g,m,T,wall,j)
        a<xj<b && throw(GeometryError("the removed open interval contains the genus-two extremum at tile $j"))
    end
    # The bisection closure checks that the specified extra gap is indeed
    # reconstructed on this wall (including its direction).
    for x in (a,b)
        h = _paper_height(g,m,T,x,wall)/T
        xhat = _paper_cut(g,m,T,wall,h)
        abs(x-xhat) <= max(1e-9*(last(g.e)-first(g.e)),8eps(maximum(abs,g.e))) ||
            throw(GeometryError("could not reconstruct the additional gap endpoint"))
    end
end

function _paper_accept(E,merged,wall,n,m,g)
    fine = _paper_geometry(g.e,2g.order)
    T = sum(g.H)/n
    T>0 && isfinite(T) || throw(GeometryError("invalid elliptic period"))
    # merged[4] or merged[5] may be trimmed by at most one tile. All other
    # branch endpoints are the supplied working-band endpoints.
    trim_below = merged[4]!=g.e[4]
    trim_above = merged[5]!=g.e[5]
    !(trim_below && trim_above) && merged[4]<=g.e[4]<g.e[5]<=merged[5] ||
        throw(GeometryError("inadmissible hidden genus-two mouth"))
    for gg in (g,fine)
        maximum(abs.([gg.left;gg.H]/T - [n,m,n-m]))<=1e-7 ||
            throw(GeometryError("genus-two heights are not quantized at the requested degree"))
        if trim_below || trim_above
            w = trim_below ? 2 : 3
            x = merged[trim_below ? 4 : 5]
            h = _paper_height(gg,m,T,x,w)/T
            amount = trim_below ? m-h : h-m
            -1e-9<=amount<=1+1e-9 || throw(GeometryError("Stiefel mouth trim exceeds one tile"))
            abs(_paper_cut(gg,m,T,w,h)-x)<=max(1e-9*(last(E)-first(E)),8eps(maximum(abs,E))) ||
                throw(GeometryError("could not reconstruct trimmed genus-two endpoint"))
        end
    end
    a,b = E[2wall:2wall+1]
    _paper_gap_check(fine,m,n,T,wall,a,b)
    omega = _paper_omega(fine)
    maximum(abs.((2/T).*omega[1,:] .- [n,n-m]))<=1e-7 ||
        throw(GeometryError("genus-two normalized periods fail quantization"))
    walls = wall==2 ? [1,2,2,3] : [1,2,3,3]
    PaperFilter(n,m,T,omega,g,copy(E),walls)
end

# Move one slit mouth into the original gap. For each bracket we check
# integral wall heights close to the observed mode before expensive inversion.
function _paper_trimmed(E,merged,wall,n,order,g0)
    q(g) = n*g.H[1]/sum(g.H)
    base = q(g0)
    lo_mode,hi_mode = max(1,floor(Int,base)-2),min(n-1,ceil(Int,base)+2)
    for index in (4,5)
        function at(s)
            e = copy(merged)
            e[index] = (1-s)*merged[index] + s*merged[9-index]
            _paper_geometry(e,order)
        end
        s0,gprev = 0.,g0
        for s1 in (.125,.25,.5,.75,.9,.975,.995,.9995,.99998)
            gnext = try
                at(s1)
            catch err
                _recoverable(err) || rethrow()
                break
            end
            x,y=q(gprev),q(gnext)
            for m in max(lo_mode,ceil(Int,min(x,y))):min(hi_mode,floor(Int,max(x,y)))
                left,right = s0,s1
                g = gnext
                for _ in 1:60
                    mid=(left+right)/2
                    g=at(mid)
                    abs(q(g)-m)<=2e-12 && break
                    if (q(g)<m)==(x<y)
                        left=mid
                    else
                        right=mid
                    end
                end
                try
                    return _paper_accept(E,merged,wall,n,m,g)
                catch err
                    _recoverable(err) || rethrow()
                end
            end
            s0,gprev = s1,gnext
        end
    end
    nothing
end

"""
    paper_fourbands(endpoints, n; order=128)

Search the article's genus-two Stiefel family for a four-band problem with
targets −,+,+,+ in the caller's real x-coordinate. Merge each adjacent pair
of positive intervals in turn. Check direct quantization, then try a single
mouth trim (7), retaining an additional open gap containing no integer-height
extremum. Return the first numerical solution with all eight working endpoints.
Failure to find this representation does not classify other genera or modes.
"""
function paper_fourbands(endpoints,n::Integer;order::Int=128)
    !(n isa Bool) && 2<=n<2^53 || throw(ArgumentError("genus-two inverse requires integer 2≤n<2^53"))
    order>=2 || throw(ArgumentError("quadrature order must be at least two"))
    E = Float64.(collect(endpoints))
    _check_endpoints(E)
    faces=Tuple{Vector{Float64},Int,Geometry}[]
    for wall in (2,3)
        merged = [E[1:2wall-1]; E[2wall+2:end]]
        try
            g=_paper_geometry(merged,order)
            push!(faces,(merged,wall,g))
            m=round(Int,n*g.H[1]/sum(g.H))
            0<m<n && return _paper_accept(E,merged,wall,Int(n),m,g)
        catch err
            _recoverable(err) || rethrow()
        end
    end
    for (merged,wall,g) in faces
        trial=_paper_trimmed(E,merged,wall,Int(n),order,g)
        isnothing(trial) || return trial
    end
    throw(GeometryError("no checked genus-two representation found by the adjacent-band merge and single-mouth trim search"))
end
