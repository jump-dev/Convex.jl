using Convex
using Base.Test

TOL = 1e-2

LPsolver() = get_default_solver()

if !can_solve_mip(get_default_solver())
	using GLPKMathProgInterface
	MIPsolver() = GLPKSolverMIP()
else
	MIPsolver() = get_default_solver()
end

@testset "Mixed Integer Programs" begin

  @testset "lp fallback interface" begin
    x = Variable()
    p = minimize(x, x>=4.3)
    @test vexity(p) == AffineVexity()
    solve!(p, LPsolver())
    @test isapprox(p.optval, 4.3, atol=TOL)

    x = Variable(2)
    p = minimize(norm(x,1), x[1]>=4.3)
    @test vexity(p) == ConvexVexity()
    solve!(p, LPsolver())
    @test isapprox(p.optval, 4.3, atol=TOL)
  end

  @testset "integer variables" begin
    x = Variable(:Int)
    p = minimize(x, x>=4.3)
    @test vexity(p) == AffineVexity()
    solve!(p, MIPsolver())
    @test isapprox(p.optval, 5, atol=TOL)

    x = Variable(2, :Int)
    p = minimize(sum(x), x>=4.3)
    @test vexity(p) == AffineVexity()
    solve!(p, MIPsolver())
    @test isapprox(p.optval, 10, atol=TOL)

    x = Variable(:Int)
    y = Variable()
    p = minimize(sum(x + y), x>=4.3, y>=7)
    @test vexity(p) == AffineVexity()
    solve!(p, MIPsolver())
    @test isapprox(p.optval, 12, atol=TOL)

    x = Variable(2, :Int)
    p = minimize(norm(x, 1), x[1]>=4.3)
    @test vexity(p) == ConvexVexity()
    solve!(p, MIPsolver())
    @test isapprox(p.optval, 5, atol=TOL)

    x = Variable(2, :Int)
    p = minimize(sum(x), x[1]>=4.3, x>=0)
    @test vexity(p) == AffineVexity()
    solve!(p, MIPsolver())
    @test isapprox(p.optval, 5, atol=TOL)

    x = Variable(2, :Int)
    p = minimize(sum(x), x>=.5)
    @test vexity(p) == AffineVexity()
    solve!(p, MIPsolver())
    @test isapprox(p.optval, 2, atol=TOL)
  end

  @testset "binary variables" begin
    x = Variable(2, :Bin)
    p = minimize(sum(x), x>=.5)
    @test vexity(p) == AffineVexity()
    solve!(p, MIPsolver())
    @test isapprox(p.optval, 2, atol=TOL)

    x = Variable(2, :Bin)
    p = minimize(sum(x), x[1]>=.5, x>=0)
    @test vexity(p) == AffineVexity()
    solve!(p, MIPsolver())
    @test isapprox(p.optval, 1, atol=TOL)
  end

end