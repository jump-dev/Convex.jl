using Convex
using FactCheck

TOL = 1e-2

# TODO: uncomment vexity checks once SDP on vars/constraints changes vexity of problem

facts("SDP Atoms") do
  context("sdp variables") do
    y = Variable((2,2), :Semidefinite)
    p = minimize(y[1,1])
    # @fact vexity(p) --> ConvexVexity()
    solve!(p)
    @fact p.optval --> roughly(0, TOL)

    y = Variable((3,3), :Semidefinite)
    p = minimize(y[1,1], y[2,2]==1)
    # @fact vexity(p) --> ConvexVexity()
    solve!(p)
    @fact p.optval --> roughly(0, TOL)

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
    @fact p.optval --> roughly(1, TOL)

    y = Variable((3, 3), :Semidefinite)
    p = minimize(trace(y), y[2,1]<=4, y[2,2]>=3)
    # @fact vexity(p) --> ConvexVexity()
    solve!(p)
    @fact p.optval --> roughly(3, TOL)

    x = Variable(Positive())
    y = Semidefinite(3)
    p = minimize(y[1, 2], y[2, 1] == 1)
    # @fact vexity(p) --> ConvexVexity()
    solve!(p)
    @fact p.optval --> roughly(1, TOL)
  end

  context("sdp constraints") do
    # This test fails on Mosek
    x = Variable(Positive())
    y = Variable((3, 3))
    p = minimize(x + y[1, 1], isposdef(y), x >= 1, y[2, 1] == 1)
    # @fact vexity(p) --> ConvexVexity()
    solve!(p)
    @fact p.optval --> roughly(1, TOL)
  end

  context("nuclear norm atom") do
    y = Semidefinite(3)
    p = minimize(nuclearnorm(y), y[2,1]<=4, y[2,2]>=3, y[3,3]<=2)
    @fact vexity(p) --> ConvexVexity()
    solve!(p)
    @fact p.optval --> roughly(3, TOL)
    @fact evaluate(nuclearnorm(y)) --> roughly(3, TOL)
  end

  context("operator norm atom") do
    y = Variable((3,3))
    p = minimize(operatornorm(y), y[2,1]<=4, y[2,2]>=3, sum(y)>=12)
    @fact vexity(p) --> ConvexVexity()
    solve!(p)
    @fact p.optval --> roughly(4, TOL)
    @fact evaluate(operatornorm(y)) --> roughly(4, TOL)
  end

  context("sigma max atom") do
    y = Variable((3,3))
    p = minimize(sigmamax(y), y[2,1]<=4, y[2,2]>=3, sum(y)>=12)
    @fact vexity(p) --> ConvexVexity()
    solve!(p)
    @fact p.optval --> roughly(4, TOL)
    @fact evaluate(sigmamax(y)) --> roughly(4, TOL)
  end

  context("lambda max atom") do
    y = Semidefinite(3)
    p = minimize(lambdamax(y), y[1,1]>=4)
    @fact vexity(p) --> ConvexVexity()
    solve!(p)
    @fact p.optval --> roughly(4, TOL)
    @fact evaluate(lambdamax(y)) --> roughly(4, TOL)
  end

  context("lambda min atom") do
    y = Semidefinite(3)
    p = maximize(lambdamin(y), trace(y)<=6)
    @fact vexity(p) --> ConvexVexity()
    solve!(p)
    @fact p.optval --> roughly(2, TOL)
    @fact evaluate(lambdamin(y)) --> roughly(2, TOL)
  end

  context("matrix frac atom") do
    x = [1, 2, 3]
    P = Variable(3, 3)
    p = minimize(matrixfrac(x, P), P <= 2*eye(3), P >= 0.5 * eye(3))
    @fact vexity(p) --> ConvexVexity()
    solve!(p)
    @fact p.optval --> roughly(7, TOL)
    @fact evaluate(matrixfrac(x, P))[1] --> roughly(7, TOL)
  end

  context("matrix frac atom both arguments variable") do
    x = Variable(3)
    P = Variable(3, 3)
    p = minimize(matrixfrac(x, P), lambdamax(P) <= 2, x[1] >= 1)
    @fact vexity(p) --> ConvexVexity()
    solve!(p)
    @fact p.optval --> roughly(.5, TOL)
    @fact evaluate(matrixfrac(x, P))[1] --> roughly(.5, TOL)
  end

  context("sum largest eigs") do
    x = Semidefinite(3)
    p = minimize(sumlargesteigs(x, 2), x >= 1)
    solve!(p)
    @fact p.optval --> roughly(3, TOL)
    @fact evaluate(x) --> roughly(ones(3, 3), TOL)

    x = Semidefinite(3)
    p = minimize(sumlargesteigs(x, 2), [x[i,:] >= i for i=1:3]...)
    solve!(p)
    @fact p.optval --> roughly(8.4853, TOL)

    x1 = Semidefinite(3)
    p1 = minimize(lambdamax(x1), x1[1,1]>=4)
    solve!(p1)

    x2 = Semidefinite(3)
    p2 = minimize(sumlargesteigs(x2, 1), x2[1,1]>=4)
    solve!(p2)

    @fact p1.optval --> roughly(p2.optval, TOL)

    x1 = Semidefinite(3)
    p1 = minimize(lambdamax(x1), [x1[i,:] >= i for i=1:3]...)
    solve!(p1)

    x2 = Semidefinite(3)
    p2 = minimize(sumlargesteigs(x2, 1), [x2[i,:] >= i for i=1:3]...)
    solve!(p2)

    @fact p1.optval --> roughly(p2.optval, TOL)

    println(p1.optval)
  end

  context("kron atom") do
    id = eye(4)
    X = Semidefinite(4)
    W = kron(id, X)
    p = maximize(trace(W), trace(X) ≤ 1)
    @fact vexity(p) --> AffineVexity()
    solve!(p)
    @fact p.optval --> roughly(4, TOL)
  end

  context("Partial trace") do
    A = Semidefinite(2)
    B = [1 0; 0 0]
    ρ = kron(B, A)
    constraints = [partialtrace(ρ, 1, [2; 2]) == [0.09942819 0.29923607; 0.29923607 0.90057181], ρ in :SDP]
    p = satisfy(constraints)
    solve!(p)
    @fact evaluate(ρ) --> roughly([0.09942819 0.29923607 0 0; 0.299237 0.900572 0 0;0 0 0 0;0 0 0 0], TOL)
    @fact evaluate(partialtrace(ρ, 1, [2; 2])) -->roughly([1.0 0;0 0], TOL)
  end
end
