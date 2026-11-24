@testset "Two seam phases with oscillations on all three bands" begin
    for n in (16, 17, 32)
        solution = oscillating_threebands(n, 0.1)
        f, E = solution.filter, solution.ends
        m = f.m
        @test f.modulus == maximum(solution.moduli)
        @test issorted(E) && E[4] < E[5]
        n == 16 && @test E[3:5] ≈ [5.9056968832e-10, 7.7036775061e-6, 1.0701086878e-5] rtol=1e-8
        # Evaluate in physical x at every interior extremum, including the
        # first broad contact at the seam. The chart values provide a second
        # calculation, independent of the inversion in welded_value.
        groups = ([(j < m ? 1 : 2, -1, j) for j in (n-1):-1:1],
                  [(1, 1, j) for j in 1:m-1], [(2, 1, j) for j in m:n-1])
        @test length.(groups) == (n - 1, m - 1, n - m)
        for (band, contacts) in enumerate(groups)
            points = [StiefelFourBands.welded_wall(f, args...) for args in contacts]
            @test all(p -> E[2band-1] < p.x < E[2band], points)
            @test issorted([p.x for p in points])
            rows = [welded_value(solution, p.x) for p in points]
            @test maximum(abs(r.value - p.value) for (r,p) in zip(rows,points)) < 1e-9
            errors = [r.normalized for r in rows]
            @test all(x -> abs(abs(x)-1) < 1e-9, errors)
            @test all(<(0), errors[2:end] .* errors[1:end-1])
        end
        @test_throws DomainError welded_value(solution, 0.0)
    end
end

@testset "Independent band weights and exchange of the two passbands" begin
    n, epsilon = 16, 0.1
    w = (2.0, 3.0, 1.0)
    direct = oscillating_threebands(n, epsilon; weights=w)
    swapped = oscillating_threebands(n, epsilon; weights=(w[1], w[3], w[2]))
    scaled = oscillating_threebands(n, 2epsilon; weights=2 .* w)
    @test scaled.ends == direct.ends
    @test !direct.reflected && swapped.reflected
    q = direct.target.q
    @test swapped.ends[4:5] ≈ q ./ reverse(direct.ends[4:5])
    for (band, chart, wall, tiles) in ((1, 1, -1, (1, 2.4, 3)),
                                      (2, 1,  1, (1, 2.4, 3)),
                                      (3, 2,  1, (9, 10.4, 11)))
        for t in tiles
            contact = StiefelFourBands.welded_wall(direct.filter, chart, wall, t)
            u = welded_value(direct, contact.x)
            v = welded_value(swapped, q/contact.x)
            @test u.value ≈ contact.value atol=1e-10
            @test v.value ≈ u.value atol=1e-10
            @test v.normalized ≈ u.normalized atol=1e-9
            @test u.normalized ≈ w[band]*(u.value-wall)/epsilon atol=1e-12
            isinteger(t) && @test abs(u.normalized) ≈ 1 atol=1e-9
        end
    end
    # The left edge of the broad band lies on the exposed seam in this example.
    x = direct.ends[5]
    @test welded_value(swapped, q/x).normalized ≈ welded_value(direct, x).normalized atol=1e-9
end

@testset "Three bands with exponentially small physical endpoints" begin
    S = StiefelFourBands
    for n in (1024, 1_000_000)
        s = oscillating_threebands(n, 0.1)
        f, E, outer = s.filter, s.ends, s.target
        @test outer.q > 0 && Float64(outer.q) == 0
        @test all(>(0), diff(E))
        @test length(f.coeff) == 28
        # Independent elliptic-integral inversion in physical x, including
        # points within a tile of either end and points near the middle.
        for v in (0.125, outer.M/3, outer.M/2, 3outer.M/4, outer.M-0.125)
            x = S._weld_abscissa(outer, v)
            k = outer.q
            d = (k/x)^2
            c2 = d * (1-x^2) / (1-k^2)
            amplitude = sqrt((1-d)/(1-k^2))
            reference = amplitude * S._zolo_rf(c2, d, one(k)) / (2outer.K)
            @test v ≈ reference atol=2e-9 rtol=0
            @test S._weld_height(outer, x) ≈ v atol=2e-9 rtol=0
        end
        for (block, wall, t) in ((1,-1,0.25), (1,1,0.75), (1,1,f.m-1.25),
                                  (2,-1,f.m+0.25), (2,1,f.m+0.5), (2,1,n-0.4))
            p = S.welded_wall(f, block, wall, t, outer)
            @test welded_value(s, p.x).value ≈ p.value atol=1e-9 rtol=0
        end
        @test all(x -> abs(welded_value(s, x).normalized) <= 1+1e-9, E)
        @test_throws DomainError welded_value(s, (E[4]+E[5])/2)
    end
end
