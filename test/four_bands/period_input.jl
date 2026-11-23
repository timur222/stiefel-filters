using LinearAlgebra: isposdef, Symmetric

@testset "Marked period input, quantization and locus" begin
    demo = load_filter(joinpath(@__DIR__,"../../data/four_bands","filter.toml"))
    p = period_parameters(demo)
    @test 2p.omega[1,:] ≈ imag(p.tau).*[p.n,p.n-p.m1,p.n-p.m2] atol=1e-10
    full = filter_from_periods(p.tau,im*p.omega; trims=p.trims,sides=p.sides)
    @test full.geom.e ≈ demo.geom.e atol=1e-10
    @test full.endpoints ≈ demo.endpoints atol=1e-10
    @test evaluate(full,-.5).direct ≈ evaluate(demo,-.5).direct atol=1e-10
    free = filter_from_periods(p.tau,p.n,p.m1,p.m2;omega22=p.omega22,omega33=p.omega33)
    @test free.geom.e ≈ demo.geom.e atol=1e-10
    @test free.diagnostics.omega23 ≈ p.omega23 atol=1e-10
    @test free.trims == (0.,0.)
    @test_throws ArgumentError filter_from_periods(1.01*p.tau,p.omega)
    @test_throws ArgumentError filter_from_periods(p.tau,p.omega + im*p.omega)
    @test_throws ArgumentError filter_from_periods(p.tau,zeros(3,3))
    @test_throws ArgumentError filter_from_periods(p.tau,ones(2,2))
    wrong=copy(p.omega)
    wrong[2,3] += .02
    wrong[3,2] += .02
    @test isposdef(Symmetric(wrong))
    @test_throws StiefelFourBands.GeometryError filter_from_periods(p.tau,wrong)
end
