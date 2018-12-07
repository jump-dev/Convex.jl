@testset "SDP and Exp Atoms: $solver" for solver in solvers
    if can_solve_sdp(solver) && can_solve_exp(solver)
        tol = 1e-2
        @testset "log det atom" begin
            x = Variable(2, 2)
            p = maximize(logdet(x), [x[1, 1] == 1, x[2, 2] == 1])
            @test vexity(p) == ConvexVexity()
            solve!(p, solver)
            @test p.optval ≈ 0 atol=tol
            @test evaluate(logdet(x)) ≈ 0 atol=tol
        end
    end
end
