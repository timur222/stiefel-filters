@testset "All-endpoint degree-three solution" begin
    for j in 1:11
        # Reproducible ordered inputs with all seven spacings changing.
        e = cumsum([0.08+3sin(j+i^2+.27)^2 for i in 1:8])
        f = cubic_fourbands(e)
        @test f.n==3 && 0<f.mu<1
        @test elliptic_constants(f.T).mu ≈ f.mu atol=2e-14
        for b in 1:4
            xs = range(e[2b-1],e[2b],length=27)
            rows = [wall_values(f,x,b) for x in xs]
            @test (rows[1].ripple,rows[end].ripple)==(-1.,1.)
            @test (rows[1].reciprocal_ripple,rows[end].reciprocal_ripple)==(1.,-1.)
            @test minimum(diff([v.direct for v in rows]))>0
            @test all(v->abs(v.ripple)<=1 && abs(v.reciprocal_ripple)<=1,rows)
            @test all(v->isapprox(v.reciprocal*v.direct,1-f.mu^2;atol=1e-14),rows)
        end
        # No normalization of r at infinity is imposed. Invariance under a
        # non-affine real coordinate change is checked off the real bands.
        m(x)=(.9x+2)/(.02x+3)
        z=.4+.3im
        moved=cubic_fourbands(m.(e))
        @test evaluate(moved,m(z)).direct ≈ evaluate(f,z).direct rtol=1e-12
    end
    # Compare with an independently prescribed cubic constructed as a
    # symmetric-arrow-matrix eigenvalue problem in the forward Pick family.
    pick = StiefelFourBands.pick_fourbands(1.02im;poles=(-.8,1.7),residues=(.7,1.2),offset=.4)
    c = cubic_fourbands(pick.endpoints)
    @test c.T ≈ pick.T atol=2e-12
    for z in (-1.7,-.01,.45,2.31,.4+.6im)
        @test evaluate(c,z).direct ≈ evaluate(pick,z).direct atol=2e-12
    end
    # μ² underflows Float64 and the mid-stopband still has a resolved ripple.
    tiny = cubic_fourbands([0.,1e-300,1e100,2e100,3e100,4e100,5e100,6e100])
    @test 0<tiny.mu<1e-180
    @test abs(wall_values(tiny,5e-301,1).ripple)<1e-14
    # L is positive even when μ rounds to one in Float64.
    large = cubic_fourbands([-1.,0.,1e-200,.3,.5,.7,.9,1.])
    @test large.mu==1.
    @test 0<wall_values(large,large.endpoints[3],2).direct<1e-198
    @test_throws ArgumentError cubic_fourbands(1:7)
    @test_throws ArgumentError cubic_fourbands(ones(8))
end
