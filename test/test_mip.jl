@testset "Mixed Integer Programs: $solver" for solver in solvers
    if can_solve_mip(solver)
        mip_solver = !can_solve_mip(solver) ? GLPKSolverMIP() : solver

        @testset "lp fallback interface" begin
            x = Variable()
            p = minimize(x, x>=4.3)
            @test vexity(p) == AffineVexity()
            solve!(p, solver)
            @test p.optval ≈ 4.3 atol=TOL

            x = Variable(2)
            p = minimize(norm(x,1), x[1]>=4.3)
            @test vexity(p) == ConvexVexity()
            solve!(p, solver)
            @test p.optval ≈ 4.3 atol=TOL
        end

        @testset "integer variables" begin
            y = Variable()
            vartype!(y, IntVar)

            for x in [ Variable(:Int), Variable(vartype = IntVar), y ]
                @test vartype(x) == IntVar
                p = minimize(x, x>=4.3)
                @test vexity(p) == AffineVexity()
                solve!(p, mip_solver)
                @test p.optval ≈ 5 atol=TOL
            end

            x = Variable(2, :Int)
            p = minimize(sum(x), x>=4.3)
            @test vexity(p) == AffineVexity()
            solve!(p, mip_solver)
            @test p.optval ≈ 10 atol=TOL

            x = Variable(:Int)
            y = Variable()
            p = minimize(sum(x + y), x>=4.3, y>=7)
            @test vexity(p) == AffineVexity()
            solve!(p, mip_solver)
            @test p.optval ≈ 12 atol=TOL

            x = Variable(2, :Int)
            p = minimize(norm(x, 1), x[1]>=4.3)
            @test vexity(p) == ConvexVexity()
            solve!(p, mip_solver)
            @test p.optval ≈ 5 atol=TOL

            x = Variable(2, :Int)
            p = minimize(sum(x), x[1]>=4.3, x>=0)
            @test vexity(p) == AffineVexity()
            solve!(p, mip_solver)
            @test p.optval ≈ 5 atol=TOL

            x = Variable(2, :Int)
            p = minimize(sum(x), x>=.5)
            @test vexity(p) == AffineVexity()
            solve!(p, mip_solver)
            @test p.optval ≈ 2 atol=TOL
        end

        @testset "binary variables" begin
            x = Variable(2, :Bin)
            p = minimize(sum(x), x>=.5)
            @test vexity(p) == AffineVexity()
            solve!(p, mip_solver)
            @test p.optval ≈ 2 atol=TOL

            x = Variable(2, :Bin)
            p = minimize(sum(x), x[1]>=.5, x>=0)
            @test vexity(p) == AffineVexity()
            solve!(p, mip_solver)
            @test p.optval ≈ 1 atol=TOL
        end
    end
end
