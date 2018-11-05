using Convex
using Test
import LinearAlgebra.I

TOL = 1e-2
eye(n) = Matrix(1.0I, n, n)


# TODO: uncomment vexity checks once SDP on vars/constraints changes vexity of problem

@testset "SDP Atoms" begin
  @testset "sdp variables" begin
    y = Variable((2,2), :Semidefinite)
    p = minimize(y[1,1])
    # @fact vexity(p) --> ConvexVexity()
    solve!(p)
    @test p.optval ≈ 0 atol=TOL

    y = Variable((3,3), :Semidefinite)
    p = minimize(y[1,1], y[2,2]==1)
    # @fact vexity(p) --> ConvexVexity()
    solve!(p)
    @test p.optval ≈ 0 atol=TOL

    # Solution is obtained as y[2,2] -> infinity
    # This test fails on Mosek. See
    # https://github.com/JuliaOpt/Mosek.jl/issues/29
    # y = Variable((2, 2), :Semidefinite)
    # p = minimize(y[1, 1], y[1, 2] == 1)
    # # @fact vexity(p) --> ConvexVexity()
    # solve!(p)
    # @fact p.optval --> roughly(0, TOL)

    y = Semidefinite(3)
    p = minimize(sum(diag(y)), y[1, 1] == 1)
    # @fact vexity(p) --> ConvexVexity()
    solve!(p)
    @test p.optval ≈ 1 atol=TOL

    y = Variable((3, 3), :Semidefinite)
    p = minimize(tr(y), y[2,1]<=4, y[2,2]>=3)
    # @fact vexity(p) --> ConvexVexity()
    solve!(p)
    @test p.optval ≈ 3 atol=TOL

    x = Variable(Positive())
    y = Semidefinite(3)
    p = minimize(y[1, 2], y[2, 1] == 1)
    # @fact vexity(p) --> ConvexVexity()
    solve!(p)
    @test p.optval ≈ 1 atol=TOL
  end

  @testset "sdp constraints" begin
    # This test fails on Mosek
    x = Variable(Positive())
    y = Variable((3, 3))
    p = minimize(x + y[1, 1], isposdef(y), x >= 1, y[2, 1] == 1)
    # @fact vexity(p) --> ConvexVexity()
    solve!(p)
    @test p.optval ≈ 1 atol=TOL
  end

  @testset "nuclear norm atom" begin
    y = Semidefinite(3)
    p = minimize(nuclearnorm(y), y[2,1]<=4, y[2,2]>=3, y[3,3]<=2)
    @test vexity(p) == ConvexVexity()
    solve!(p)
    @test p.optval ≈ 3 atol=TOL
    @test evaluate(nuclearnorm(y)) ≈ 3 atol=TOL
  end

  @testset "operator norm atom" begin
    y = Variable((3,3))
    p = minimize(operatornorm(y), y[2,1]<=4, y[2,2]>=3, sum(y)>=12)
    @test vexity(p) == ConvexVexity()
    solve!(p)
    @test p.optval ≈ 4 atol=TOL
    @test evaluate(operatornorm(y)) ≈ 4 atol=TOL
  end

  @testset "sigma max atom" begin
    y = Variable((3,3))
    p = minimize(sigmamax(y), y[2,1]<=4, y[2,2]>=3, sum(y)>=12)
    @test vexity(p) == ConvexVexity()
    solve!(p)
    @test p.optval ≈ 4 atol=TOL
    @test evaluate(sigmamax(y)) ≈ 4 atol=TOL
  end

  @testset "lambda max atom" begin
    y = Semidefinite(3)
    p = minimize(lambdamax(y), y[1,1]>=4)
    @test vexity(p) == ConvexVexity()
    solve!(p)
    @test p.optval ≈ 4 atol=TOL
    @test evaluate(lambdamax(y)) ≈ 4 atol=TOL
  end

  @testset "lambda min atom" begin
    y = Semidefinite(3)
    p = maximize(lambdamin(y), tr(y)<=6)
    @test vexity(p) == ConvexVexity()
    solve!(p)
    @test p.optval ≈ 2 atol=TOL
    @test evaluate(lambdamin(y)) ≈ 2 atol=TOL
  end

  @testset "matrix frac atom" begin
    x = [1, 2, 3]
    P = Variable(3, 3)
    p = minimize(matrixfrac(x, P), P <= 2*eye(3), P >= 0.5 * eye(3))
    @test vexity(p) == ConvexVexity()
    solve!(p)
    @test p.optval ≈ 7 atol=TOL
    @test (evaluate(matrixfrac(x, P)))[1] ≈ 7 atol=TOL
  end

  @testset "matrix frac atom both arguments variable" begin
    x = Variable(3)
    P = Variable(3, 3)
    p = minimize(matrixfrac(x, P), lambdamax(P) <= 2, x[1] >= 1)
    @test vexity(p) == ConvexVexity()
    solve!(p)
    @test p.optval ≈ 0.5 atol=TOL
    @test (evaluate(matrixfrac(x, P)))[1] ≈ 0.5 atol=TOL
  end

  @testset "sum largest eigs" begin
    x = Semidefinite(3)
    p = minimize(sumlargesteigs(x, 2), x >= 1)
    solve!(p)
    @test p.optval ≈ 3 atol=TOL
    @test evaluate(x) ≈ ones(3, 3) atol=TOL

    x = Semidefinite(3)
    p = minimize(sumlargesteigs(x, 2), [x[i,:] >= i for i=1:3]...)
    solve!(p)
    @test p.optval ≈ 8.4853 atol=TOL

    x1 = Semidefinite(3)
    p1 = minimize(lambdamax(x1), x1[1,1]>=4)
    solve!(p1)

    x2 = Semidefinite(3)
    p2 = minimize(sumlargesteigs(x2, 1), x2[1,1]>=4)
    solve!(p2)

    @test p1.optval ≈ p2.optval atol=TOL

    x1 = Semidefinite(3)
    p1 = minimize(lambdamax(x1), [x1[i,:] >= i for i=1:3]...)
    solve!(p1)

    x2 = Semidefinite(3)
    p2 = minimize(sumlargesteigs(x2, 1), [x2[i,:] >= i for i=1:3]...)
    solve!(p2)

    @test p1.optval ≈ p2.optval atol=TOL

    println(p1.optval)
  end

  @testset "kron atom" begin
    id = eye(4)
    X = Semidefinite(4)
    W = kron(id, X)
    p = maximize(tr(W), tr(X) ≤ 1)
    @test vexity(p) == AffineVexity()
    solve!(p)
    @test p.optval ≈ 4 atol=TOL
  end

  @testset "Partial trace" begin
    A = Semidefinite(2)
    B = [1 0; 0 0]
    ρ = kron(B, A)
    constraints = [partialtrace(ρ, 1, [2; 2]) == [0.09942819 0.29923607; 0.29923607 0.90057181], ρ in :SDP]
    p = satisfy(constraints)
    solve!(p)
    @test evaluate(ρ) ≈ [0.09942819 0.29923607 0 0; 0.299237 0.900572 0 0; 0 0 0 0; 0 0 0 0] atol=TOL
    @test evaluate(partialtrace(ρ, 1, [2; 2])) ≈ [1.0 0; 0 0] atol=TOL
  end
end
