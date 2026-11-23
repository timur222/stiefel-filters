# Reference periods and U values use Python SC integration on saved
# prevertices from the degree-three rounded working-endpoint example.
@testset "Genus-three theta periods, shifts, and inverse" begin
    e = [-1., 0., .10525, .10772, .3515875750036856, .37827, .9341125184083909, 1.]
    curve = StiefelFourBands.ThetaCurve(StiefelFourBands.geometry(e,128))
    reference = [1.0500877568409939 .700058504560662 .3500292522803307;
                 .700058504560662 .8522154941218153 .3948353047488519;
                 .3500292522803307 .3948353047488519 .5744824994896568]
    @test maximum(abs, curve.omega-reference) < 5e-12
    @test StiefelFourBands.riemann_constant_residual(curve) < 1e-10
    z = [.2+.1im, -.17+.2im, .1-.14im]
    val, grad = StiefelFourBands.theta_value_gradient(curve, z, ([1,0,1], [0,1,1]))
    translated, _ = StiefelFourBands.theta_value_gradient(curve, z .+ [1,0,0], ([1,0,1], [0,1,1]))
    @test translated ≈ -val atol=5e-12
    expanded, _ = StiefelFourBands.theta_value_gradient(curve, z, ([101,0,1], [2,1,1]))
    @test expanded ≈ -val atol=5e-12
    for (point, u) in [
        (.3+.2im, .15813840914874178+1.1201082169195882im),
        (-.4+.35im, -.6439571429351623+1.4629067494109718im),
        (.8+.15im, .4291634351696122+1.5379488802551262im)]
        @test StiefelFourBands.theta_coordinate(curve, point) ≈ u atol=1e-8
    end
end

@testset "Controlled genus-three theta truncation" begin
    e = [-1., 0., .10525, .10772, .3515875750036856, .37827, .9341125184083909, 1.]
    c = StiefelFourBands.ThetaCurve(StiefelFourBands.geometry(e,128))
    z = [.28+.71im, -.37+.2im, .41-.6im]
    ch = ([101,0,1], [2,1,1])
    controlled = StiefelFourBands.theta_controlled(c,z,ch; tolerance=1e-11)
    @test controlled.lambda_lower > 0
    @test max(controlled.log_value_tail,controlled.log_gradient_tail) <= log(1e-11)
    enlarged = StiefelFourBands.ThetaCurve(StiefelFourBands.geometry(e,128); cutoff=controlled.cutoff+2)
    v, grad = StiefelFourBands.theta_value_gradient(enlarged,z,ch; scaled=true)
    @test abs(controlled.value-v) < 1e-12
    @test maximum(abs,controlled.gradient-grad) < 1e-11
end
