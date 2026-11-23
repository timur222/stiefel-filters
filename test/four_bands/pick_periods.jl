@testset "Marked genus-five periods of a cubic Pick pullback" begin
    S=StiefelFourBands
    f=pick_fourbands(.8im;poles=(-1.1,1.7),residues=(.7,1.3),offset=.2)
    Ω=pick_period_matrix(f;order=96)
    @test size(Ω)==(5,5) && S.isposdef(S.Symmetric(Ω))
    @test Ω ≈ pick_period_matrix(f;order=192) atol=2e-11
    # The even null K+[2] vanishes on this marked hyperelliptic Jacobian:
    # unlike an odd characteristic, it is not zero on every period matrix.
    @test abs(θ[2,3,5,7,9,11](im*Ω)) < 4e-13
    @test Ω*[1.,-1.,1.,-1.,1.] ≈ .4*[1.,0.,1.,0.,1.] atol=2e-12
    g=pick_from_periods(.8im,im*Ω;order=96)
    @test collect(g.poles) ≈ [-1.4,1.4] atol=3e-6
    @test collect(g.residues) ≈ [.7,1.3] atol=3e-6
    @test g.offset ≈ .5 atol=3e-6
    @test (g.endpoints .+ .3) ≈ f.endpoints atol=3e-6
    @test pick_period_matrix(g) ≈ Ω atol=4e-9
    @test evaluate(g,.4+.2im).direct ≈ evaluate(f,.7+.2im).direct atol=3e-6
    # Positive definite perturbation in the complement of the Pick locus.
    v=[1.,1.,0.,0.,0.]
    bad=Ω+.02*v*transpose(v)
    @test_throws S.GeometryError pick_from_periods(.8im,bad;order=64)
    @test_throws ArgumentError pick_from_periods(.8im,zeros(5,5))
    @test_throws ArgumentError pick_from_periods(.6im,Ω)
end

@testset "Genus-five indexed characteristics and theta sum" begin
    @test branch_characteristic(5,(3,5,7,9,11)) == ([1,1,1,1,1],[1,0,1,0,1])
    @test branch_characteristic(5,(12,)) == ([0,0,0,0,0],[1,1,1,1,1])
    @test_throws ArgumentError branch_characteristic(5,(13,))
    # A diagonal Π factors exactly into five ordinary elliptic theta series.
    diagvals=[.8,1.1,1.4,1.,1.3]
    Π=im*StiefelFourBands.diagm(diagvals)
    z=ComplexF64[.1+.3im,-.23,.04-.1im,.15+.08im,.29]
    a,b=branch_characteristic(5,(3,7,11))
    ref=prod(sum(exp(-π*t*(k+aa/2)^2+2π*im*(k+aa/2)*(v+bb/2)) for k in -12:12)
             for (t,v,aa,bb) in zip(diagvals,z,a,b))
    @test θ[3,7,11](z,Π) ≈ ref atol=2e-14
    # Off-diagonal couplings and reduction: compare to a rectangular sum.
    Π[1,2]=Π[2,1]=.2+.1im
    @test θ[3,7,11](z,Π) ≈ theta_digits(Π,(3,7,11),z;cutoff=4) atol=3e-14
end

@testset "Branch-period marking across genera" begin
    S=StiefelFourBands
    f=paper_filter(.7im,3,1,1.25)
    @test S._branch_omega(f.geom.e,96) ≈ f.omega atol=3e-12
    e=[-1.,0.,.10525,.10772,.35159,.37827,.93411,1.]
    @test S._branch_omega(e,96) ≈ ThetaCurve(S.geometry(e,96)).omega atol=3e-12
    @test S._branch_omega([-1.,-.5,.5,1.],96)[1] > 0
end
