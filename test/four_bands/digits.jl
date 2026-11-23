@testset "Paper branch digit notation and genus-two periods" begin
    table = (([0,0],[0,0]),([1,0],[0,0]),([1,0],[1,0]),
             ([0,1],[1,0]),([0,1],[1,1]),([0,0],[1,1]))
    for i in 1:6
        @test branch_characteristic(2,i) == table[i]
    end
    @test branch_characteristic(2,"[35]") == ([1,1],[0,1])
    @test branch_characteristic(2,"3355") == ([0,0],[0,0])
    @test branch_characteristic(3,357) == (StiefelFourBands.RIEMANN_A, StiefelFourBands.RIEMANN_B)
    for i in 1:8
        @test branch_characteristic(3,i) ==
            (StiefelFourBands.BRANCH_ALPHA[i,:], StiefelFourBands.BRANCH_BETA[i,:])
        @test branch_characteristic(3,(3,5,7,i)) == StiefelFourBands.theta_char(i)
    end
    Π = paper_period_matrix(.7im,3,1,1.25)
    @test Π ≈ im*[1.05 .7; .7 1.25]
    @test abs(θ[35](Π)) < 2e-15
    z = [.23+.17im, -.37+.1im]
    @test θ[35](z,Π) ≈ θ[3,5](z,Π)
    @test θ[35](-z,Π) ≈ -θ[35](z,Π)
    # Unreduced definition (19), independent of the matrix/reindexing code.
    general = Π + [.21 .09; .09 -.17]
    a,b = table[3]
    zt = z + general*[1,-1] + [2,1]
    reference = sum(begin
        v = [p,q] + a/2
        exp(π*im*(sum(v .* (general*v))+2sum(v .* (zt+b/2))))
    end for p in -13:13, q in -13:13)
    @test θ[3](zt,general) ≈ reference rtol=1e-12
    x = paper_branchpoints(Π)
    @test x[[1,2,6]] == [0.,1.,Inf]
    @test 1 < x[3] < x[4] < x[5]
    for s in (2,3,4,5)
        a,b = table[s]
        u = (Π*a+b)/2
        @test abs(θ[35](u,Π)) < 1e-13
        # (21) has a removable 0/0 at p₄. Switch there to (20) with k=3.
        ratio = s==4 ? θ[625](Π)*θ[15](u,Π)/(θ[125](Π)*θ[65](u,Π)) :
                         θ[1](Π)*θ[26](u,Π)/(θ[6](Π)*θ[2](u,Π))
        @test -ratio^2 ≈ x[s] rtol=2e-13
    end
    @test_throws ArgumentError theta_digits(Π,"7",z)
    @test_throws ArgumentError theta_digits(Π,"3a",z)
    @test_throws ArgumentError theta_digits(im*[1. 1.; 1. 1.],1,z)
    @test_throws ArgumentError paper_period_matrix(.7,3,1,1.3)
    @test_throws ArgumentError paper_period_matrix(.7im,3,1,.6)
    @test_throws ArgumentError paper_branchpoints(im*[1. -.2; -.2 2.])
end

@testset "Genus-two SC coordinate against the paper's theta projection" begin
    S = StiefelFourBands
    f = paper_filter(.7im,3,1,1.25)
    e, M = f.geom.e, f.geom.mom
    Π = im*f.omega
    A = 2M[[2,4],1:2]
    A[2,:] .*= -1
    z = .35+.2im
    d = z-e[3]
    abelian = [S._rule(0.,1.,96) do t
        x = e[3]+d*t^2
        2sqrt(d)*x^k/prod(sqrt(x-e[j]) for j in (1,2,4,5,6))
    end for k in 0:1]
    a,b = branch_characteristic(2,3)
    w0 = (b-Π*a)/2
    w = w0 + (transpose(A) \ abelian)
    @test abs(θ[35](w,Π)) < 5e-14
    @test -(θ[1](Π)*θ[26](w,Π)/(θ[6](Π)*θ[2](w,Π)))^2 ≈ (1+z)/(1-z) atol=2e-13
    @test 1+4(w[1]-w0[1]) ≈ S._upper_coordinate(f.geom,z) atol=5e-13
end
