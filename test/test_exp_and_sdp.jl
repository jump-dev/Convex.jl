using Convex
using Test

TOL = 1e-2

@testset "SDP and Exp Atoms" begin

  @testset "log det atom" begin
    x = Variable(2, 2)
    p = maximize(logdet(x), [x[1, 1] == 1, x[2, 2] == 1])
    @test vexity(p) == ConvexVexity()
    solve!(p)
    @test isapprox(p.optval, 0, atol=TOL)
    @test isapprox(evaluate(logdet(x)), 0, atol=TOL)
  end

end
