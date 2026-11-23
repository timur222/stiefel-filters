include("../../src/four_bands/StiefelFourBands.jl")
using .StiefelFourBands
using Test
using TOML

include("evaluate.jl")
include("theta.jl")
include("geometry.jl")

@testset "Reusable TOML and offline plot" begin
    E=[-1.,0.,.10525,.10772,.35397,.37827,.94801,1.]
    m=synthesize(E,3)
    @test (m.m1,m.m2)==(1,2)
    @test elliptic_constants(m.T).mu ≈ .04496701750264 atol=2e-10
    mktempdir() do dir
        p=save_filter(joinpath(dir,"demo.toml"),m)
        loaded=load_filter(p)
        @test loaded.geom.e == m.geom.e
        @test evaluate(loaded,-.5).direct ≈ -.9578378141771894 atol=1e-9
        @test evaluate(loaded,.3+.2im).reciprocal ≈ .8080660716237286+.9963257216072002im atol=2e-9
        html=read(plot_filter(joinpath(dir,"demo.html"),loaded;samples=32),String)
        @test occursin("Plotly.newPlot",html)
        @test !occursin(r"<script\s+src=",html)
        spec=StiefelFourBands._figure(m;samples=32)
        @test length(spec["data"])==22
        @test all(tr->all(v->v===nothing || isfinite(v),tr["y"]) && issorted(tr["x"]),spec["data"])
        d=TOML.parsefile(p); d["T"]*=1.1
        open(io->TOML.print(io,d),p,"w")
        @test_throws Exception load_filter(p)
    end
end

include("digits.jl")
include("period_input.jl")
include("other_bands.jl")
include("cubic_inverse.jl")
include("pick_periods.jl")
include("paper_filter.jl")
include("paper_inverse.jl")

@testset "Period input files and independently specified full matrix" begin
    example(name)=joinpath(@__DIR__,"../../data/four_bands",name)
    m=read_period_filter(example("periods-g3.toml"))
    @test m.endpoints ≈ load_filter(example("filter.toml")).endpoints atol=2e-11
    @test read_period_filter(example("periods-g2.toml")) isa PaperFilter
    @test read_period_filter(example("periods-g1.toml")) isa ZolotarevFourBands
    @test source_genus(read_period_filter(example("periods-g4.toml")))==4
    @test source_genus(read_period_filter(example("periods-g5.toml")))==5
    @test source_genus(read_period_filter(example("matrix-g5.toml")))==5
    @test cubic_fourbands(-4:3).mu == 9/16
    mktempdir() do dir
        path=joinpath(dir,"input.toml")
        p=period_parameters(m)
        d=Dict("genus"=>3,"tau_im"=>imag(p.tau),"omega"=>collect(eachrow(p.omega)))
        open(io->TOML.print(io,d),path,"w")
        @test read_period_filter(path).geom.e ≈ m.geom.e atol=2e-10
        d["omega33"]=1
        open(io->TOML.print(io,d),path,"w")
        @test_throws ArgumentError read_period_filter(path)
    end
end

@testset "Complete marked matrix entry point" begin
    f=zolotarev_fourbands(.7im,3)
    Y=StiefelFourBands._branch_omega(f.endpoints[[1,2,3,8]],128)
    @test only(Y) ≈ 3*.7/2 atol=2e-12
    @test filter_from_periods(.7im,im*Y).endpoints ≈ f.endpoints
    g=filter_from_periods(.7im,im*[1.05 .7;.7 1.8];split=(2,1,.25,.65))
    @test g.endpoints ≈ read_period_filter(joinpath(@__DIR__,"../../data/four_bands","periods-g2.toml")).endpoints
    L=folded_fourbands(.8im;scale=1.2,turning=2.1,odd_root=-2.7)
    @test filter_from_periods(.8im,folded_period_matrix(L)).endpoints ≈ L.endpoints atol=1e-7
    @test read_period_filter(joinpath(@__DIR__,"../../data/four_bands","matrix-g4.toml")).endpoints ≈ L.endpoints atol=1e-7
    @test_throws ArgumentError filter_from_periods(.7im,reshape([.73],1,1))
end
