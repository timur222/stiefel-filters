@testset "Genus two from article periods and an alternation-preserving gap" begin
    f=paper_filter(.7im,3,1,1.8;split=(2,1,.25,.65))
    base=paper_filter(.7im,3,1,1.8)
    @test length(base.endpoints)==6
    @test f.endpoints[[1,2,3,4,5,8]] ≈ base.endpoints
    @test f.geom.H ≈ [.7,1.4] atol=2e-12
    @test f.geom.left ≈ 2.1 atol=2e-12
    @test f.endpoints[6:7] ≈ [.431946974933159,.552794837955880] atol=2e-12
    @test StiefelFourBands.height_tiles(f,f.endpoints[6],3) ≈ 1.25 atol=1e-12
    @test StiefelFourBands.height_tiles(f,f.endpoints[7],4) ≈ 1.65 atol=1e-12
    @test evaluate(f,.3+.2im).direct ≈ evaluate(base,.3+.2im).direct
    for x in (-2.,-.4,.08,.31,.376,.5,1.25)
        @test evaluate(f,x).direct ≈ real(evaluate(f,x+1e-7im).direct) rtol=2e-6 atol=5e-8
    end
    @test length(StiefelFourBands._figure(base;samples=16)["data"])==16
    @test length(StiefelFourBands._figure(f;samples=16)["data"])==22
    @test_throws ArgumentError paper_filter(.7im,3,1,1.8;split=(2,0,.25,.65))
end
