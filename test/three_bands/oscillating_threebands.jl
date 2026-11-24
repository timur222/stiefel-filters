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
