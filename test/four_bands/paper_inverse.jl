@testset "Genus-two inversion from eight physical endpoints" begin
    for split in ((1,0,.15,.54),(2,2,.12,.63))
        src = paper_filter(.7im,3,1,1.25;split)
        f = paper_fourbands(src.endpoints,3)
        @test f.endpoints == src.endpoints
        @test f.walls == src.walls
        @test f.omega ≈ src.omega atol=2e-10
        @test f.m==src.m && f.T≈src.T
        mapx(x)=(3x+1)/(.09x+2)
        moved = paper_fourbands(mapx.(src.endpoints),3)
        @test moved.omega ≈ f.omega atol=2e-10
        @test evaluate(moved,mapx(-.2)).direct ≈ evaluate(src,-.2).direct atol=2e-10
    end
    src = paper_filter(.7im,3,1,1.25;split=(1,0,.15,.54))
    E = copy(src.endpoints)
    E[7] = StiefelFourBands._paper_cut(src.geom,src.m,src.T,3,1.3)
    trimmed = paper_fourbands(E,3)
    @test trimmed.endpoints == E
    @test trimmed.T ≈ src.T atol=1e-9
    @test trimmed.omega ≈ src.omega atol=1e-9
    # In this representation the open gap crosses a strict interior
    # alternation point of the underlying positive wall.
    a=StiefelFourBands._paper_cut(src.geom,src.m,src.T,3,1.6)
    b=StiefelFourBands._paper_cut(src.geom,src.m,src.T,3,2.3)
    @test_throws StiefelFourBands.GeometryError StiefelFourBands._paper_gap_check(
        src.geom,src.m,src.n,src.T,3,a,b)
end
