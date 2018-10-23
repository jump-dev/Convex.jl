using Convex
using Test

TOL = 1e-3

@testset "Exp Atoms" begin

  @testset "exp atom" begin
    y = Variable()
    p = minimize(exp(y), y>=0)
    @test vexity(p) == ConvexVexity()
    solve!(p)
    @test isapprox(p.optval, 1, atol=TOL)
    @test isapprox(evaluate(exp(y)), 1, atol=TOL)

    y = Variable()
    p = minimize(exp(y), y>=1)
    @test vexity(p) == ConvexVexity()
    solve!(p)
    @test isapprox(p.optval, exp(1), atol=TOL)
    @test isapprox(evaluate(exp(y)), exp(1), atol=TOL)

    y = Variable(5)
    p = minimize(sum(exp(y)), y>=0)
    @test vexity(p) == ConvexVexity()
    solve!(p)
    @test isapprox(p.optval, 5, atol=TOL)
    @test isapprox(evaluate(sum(exp(y))), 5, atol=TOL)

    y = Variable(5)
    p = minimize(sum(exp(y)), y>=0)
    @test vexity(p) == ConvexVexity()
    solve!(p)
    @test isapprox(p.optval, 5, atol=TOL)
  end

  @testset "log atom" begin
    y = Variable()
    p = maximize(log(y), y<=1)
    @test vexity(p) == ConvexVexity()
    solve!(p)
    @test isapprox(p.optval, 0, atol=TOL)

    y = Variable()
    p = maximize(log(y), y<=2)
    @test vexity(p) == ConvexVexity()
    solve!(p)
    @test isapprox(p.optval, log(2), atol=TOL)

    y = Variable()
    p = maximize(log(y), [y<=2, exp(y)<=10])
    @test vexity(p) == ConvexVexity()
    solve!(p)
    @test isapprox(p.optval, log(2), atol=TOL)
  end

  @testset "log sum exp atom" begin
    y = Variable(5)
    p = minimize(logsumexp(y), y>=1)
    @test vexity(p) == ConvexVexity()
    solve!(p)
    @test isapprox(p.optval, log(exp(1) * 5), atol=TOL)
  end

  @testset "logistic loss atom" begin
    y = Variable(5)
    p = minimize(logisticloss(y), y>=1)
    @test vexity(p) == ConvexVexity()
    solve!(p)
    @test isapprox(p.optval, log(exp(1) + 1) * 5, atol=TOL)
  end

  @testset "entropy atom" begin
    y = Variable(5, Positive())
    p = maximize(entropy(y), sum(y)<=1)
    @test vexity(p) == ConvexVexity()
    solve!(p)
    @test isapprox(p.optval, -(log(1 / 5)), atol=TOL)
  end

  @testset "relative entropy atom" begin
    x = Variable(1)
    y = Variable(1)
    # x log (x/y)
    p = minimize(relative_entropy(x,y), y==1, x >= 2)
    @test vexity(p) == ConvexVexity()
    solve!(p)
    @test isapprox(p.optval, 2 * log(2), atol=TOL)
  end

  @testset "log perspective atom" begin
    x = Variable(1)
    y = Variable(1)
    # y log (x/y)
    p = maximize(log_perspective(x,y), y==5, x <= 10)
    @test vexity(p) == ConvexVexity()
    solve!(p)
    @test isapprox(p.optval, 5 * log(2), atol=TOL)
  end

end
