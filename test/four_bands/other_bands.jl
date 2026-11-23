@testset "Degree two elliptic family on four prescribed bands" begin
    f=zolotarev_fourbands(.7im,2)
    @test f.endpoints[1]<f.endpoints[2]<f.endpoints[3]<f.endpoints[8]
    k=f.q
    mu=elliptic_constants(f.T).mu
    # For degree two the reciprocal fraction has an elementary independent formula.
    for x in (-2.,-.8,-sqrt(k),-.2,.16,.56,1.4,.25+.37im,-.2-.34im,2+.1im)
        v=evaluate(f,x)
        exact=2*(x+k/x)/(1+sqrt(k))^2
        @test v.reciprocal ≈ exact rtol=2e-12
        @test v.direct ≈ (1-mu^2)/exact rtol=2e-12
    end
    b=f.endpoints
    @test_throws ArgumentError zolotarev_fourbands(.7im,2;gaps=((.2,1.1),(1.3,1.5)))
    @test_throws ArgumentError zolotarev_fourbands([b[1:3]...,(k+sqrt(k))/2,
        (sqrt(k)+1)/2,.9,.95,1.],2)
    # Projective covariance for non-symmetric physical intervals.
    M(x)=(3x+1)/(.2x+2)
    h=zolotarev_fourbands(M.(b),2)
    @test h.T ≈ f.T rtol=2e-13
    for x in (-.62,.2,.7,.2-.42im)
        @test evaluate(h,M(x)).direct ≈ evaluate(f,x).direct rtol=2e-12
    end
    errors=[evaluate(f,x).ripple for x in (-1.,-sqrt(k),-k,k,sqrt(k),1.)]
    @test errors ≈ [1.,-1.,1.,-1.,1.,-1.] atol=2e-12
    for i in 1:4
        vals=[wall_values(f,x,i) for x in range(b[2i-1],b[2i];length=15)]
        @test maximum(v->abs(v.ripple),vals)<=1+2e-13
    end
end

@testset "Elliptic construction at degrees 1, 5 and 3000" begin
    for (n,T) in ((1,.85),(5,.37),(3000,.0004))
        f=zolotarev_fourbands(im*T,n)
        g=zolotarev_fourbands(f.endpoints,n)
        @test g.T ≈ T rtol=2e-12
        @test evaluate(g,(g.endpoints[1]+g.endpoints[2])/2).direct ≈
              evaluate(f,(f.endpoints[1]+f.endpoints[2])/2).direct atol=1e-12
    end
end

@testset "A four-band local solution of canonical genus five" begin
    f=pick_fourbands(.7im)
    mu=elliptic_constants(f.T).mu
    @test length(unique(pick_branchpoints(f))) == 12
    @test f.n==3
    @test length(f.endpoints)==8
    for b in 1:4
        left,right=f.endpoints[2b-1:2b]
        vl,vr=wall_values(f,left,b),wall_values(f,right,b)
        @test [vl.ripple,vr.ripple] ≈ [-1.,1.] atol=3e-13
        @test [vl.reciprocal_ripple,vr.reciprocal_ripple] ≈ [1.,-1.] atol=3e-13
        @test maximum(x->abs(wall_values(f,x,b).ripple),range(left,right;length=20))<1+3e-13
    end
    for x in (-2.,.4,-.5+.4im)
        @test evaluate(f,x).direct ≈ (x^3-3x)/(x^2-1) rtol=2e-14
        @test evaluate(f,x).direct*evaluate(f,x).reciprocal ≈ 1-mu^2
    end
    # Every real v has one simple preimage on each open pole branch.
    # An independent count for the default cubic uses the discriminant of
    # x³-v x²-3x+v, equal to 4(v⁴+9v²+27).
    for v in (-1-mu,-1+mu,1-mu,1+mu)
        @test 4*(v^4+9v^2+27)>0
    end
    other=pick_fourbands(.55im;poles=(-2.,.3),residues=(.4,1.7),offset=.3)
    @test other.endpoints!=f.endpoints
    @test wall_values(other,other.endpoints[5],3).ripple ≈ -1 atol=1e-12
    @test_throws ArgumentError pick_fourbands(.7im;residues=(-1.,1.))
end

@testset "Genus four from a single double contact" begin
    for T in (.6,.95)
        f=folded_fourbands(im*T)
        @test length(unique(f.branchpoints)) == 10
        @test issorted(f.endpoints)
        e=f.endpoints
        turning=-1/(2-f.cut)
        @test e[1]<turning<e[2]
        @test e[4]<0<e[5]
        @test_throws DomainError evaluate(f,0.)
        extremals=[(e[1],1),(turning,1),(e[2],1),(e[3],2),(e[4],2),
                   (e[5],3),(e[6],3),(e[8],4)]
        values=[wall_values(f,x,b) for (x,b) in extremals]
        alt=[1.,-1.,1.,-1.,1.,-1.,1.,-1.]
        @test [v.ripple for v in values] ≈ alt atol=4e-12
        @test [v.reciprocal_ripple for v in values] ≈ -alt atol=4e-12
        for b in 1:4
            @test maximum(x->abs(wall_values(f,x,b).ripple),
                range(e[2b-1],e[2b];length=100))<=1+4e-12
        end
        z=2+.3im
        x=f.cut-1/z
        @test evaluate(f,z).direct ≈ x-9/(x+1)+2/(x-1)-2-elliptic_constants(T).mu
    end
end

@testset "Marked period input for the genus-four cubic" begin
    f=folded_fourbands(.7im;scale=1.25,turning=1.7,odd_root=-2.4)
    mu=elliptic_constants(f.T).mu
    # Recover the four raw levels in the physical marking without using
    # endpoint labels; only the lower stop level has one odd preimage.
    values=sort([evaluate(f,w).direct for w in f.branchpoints])
    @test values ≈ [-1-mu; fill(-1+mu,3); fill(1-mu,3); fill(1+mu,3)] atol=3e-13
    P=folded_period_matrix(f;order=128)
    @test size(P)==(4,4) && isposdef(Symmetric(P))
    g=folded_from_periods(.7im,im*P;order=128)
    @test g.endpoints ≈ f.endpoints atol=1e-7
    @test [g.scale,g.turning,g.odd_root] ≈ [1.25,1.7,-2.4] atol=1e-7
    @test folded_period_matrix(g;order=128) ≈ P atol=1e-9
    off=copy(P)
    off[4,4]+=.05
    @test isposdef(Symmetric(off))
    @test_throws Exception folded_from_periods(.7im,off;order=128)
    @test_throws Exception folded_from_periods(.74im,P;order=128)
end
