@testset "Elliptic values and stable wall ripple" begin
    # Independent reference values from 80-digit mpmath Jacobi sn.
    for (T, μ, sn) in (
        (0.7, .04495017213363394326, .51711711985105769+.32807897456063932im),
        (0.25, .000013949369424157398, 1.0714698988349161+.12991729751636966im),
        (0.05, 2.0631600250171361e-27, 1.0000011498192051-.00000353880094387021im))
        c = StiefelFourBands.elliptic_constants(T)
        @test c.mu ≈ μ rtol=5e-15
        @test StiefelFourBands.elliptic_value(.21+.17im,T) ≈ sn atol=1e-14
        @test StiefelFourBands.boundary_ripple(T,T,1) ≈ 1 atol=2e-14
        @test StiefelFourBands.boundary_ripple(0,T,-1) ≈ 1 atol=2e-14
    end
    c = StiefelFourBands.elliptic_constants(0.0005)
    @test iszero(c.mu)
    @test c.log_mu ≈ log(4)-π/.0005
    @test StiefelFourBands.boundary_ripple(.173*.0005,.0005) ≈ -cospi(.173)
    @test_throws DomainError StiefelFourBands.elliptic_value(.7im,.7)
end

@testset "Direct and reciprocal fractions" begin
    E = [-1., 0., .10525, .10772, .35397, .37827, .94801, 1.]
    f = StiefelFourBands.synthesize(E, 3; levels=(1,2), sides=(:above,:above))
    c = StiefelFourBands.elliptic_constants(f.T)
    for (x, ref) in ((-.5, -.9578378141771894), (.3+.2im, .49004292140356626-.604210948134688im))
        v = StiefelFourBands.evaluate(f,x)
        @test v.direct ≈ ref atol=1e-8
        @test v.reciprocal*v.direct ≈ (1-c.mu^2) atol=1e-13
        @test StiefelFourBands.evaluate(f,conj(x)).direct ≈ conj(v.direct)
    end
    for (j, (a,b)) in enumerate(zip(E[1:2:end], E[2:2:end]))
        for x in (a, (a+b)/2, b)
            v = StiefelFourBands.wall_values(f, x, j)
            @test abs(v.ripple) <= 1 + 1e-10
            @test (v.reciprocal-v.target)/c.mu ≈ v.reciprocal_ripple atol=1e-12
        end
    end
    # One interior point on each horizontal trace and both exterior rays.
    for x in (-2., .05, .23, .57, 1.2)
        boundary = StiefelFourBands.evaluate(f,x)
        upper = StiefelFourBands.evaluate(f,x+1e-7im)
        @test boundary.direct ≈ real(upper.direct) rtol=3e-7 atol=3e-8
        @test boundary.reciprocal ≈ real(upper.reciprocal) rtol=3e-7 atol=3e-8
    end
end
