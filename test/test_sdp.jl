using Convex
using FactCheck

TOL = 1e-2

# TODO: uncomment vexity checks once SDP on vars/constraints changes vexity of problem

facts("SDP Atoms") do
  context("sdp variables") do
    y = Variable((2,2), :Semidefinite)
    p = minimize(y[1,1])
    # @fact vexity(p) => ConvexVexity()
    solve!(p)
    @fact p.optval => roughly(0, TOL)

    y = Variable((3,3), :Semidefinite)
    p = minimize(y[1,1], y[2,2]==1)
    # @fact vexity(p) => ConvexVexity()
    solve!(p)
    @fact p.optval => roughly(0, TOL)

    # Solution is obtained as y[2,2] -> infinity
    # This test fails on Mosek. See
    # https://github.com/JuliaOpt/Mosek.jl/issues/29
    # y = Variable((2, 2), :Semidefinite)
    # p = minimize(y[1, 1], y[1, 2] == 1)
    # # @fact vexity(p) => ConvexVexity()
    # solve!(p)
    # @fact p.optval => roughly(0, TOL)

    y = Semidefinite(3)
    p = minimize(sum(diag(y)), y[1, 1] == 1)
    # @fact vexity(p) => ConvexVexity()
    solve!(p)
    @fact p.optval => roughly(1, TOL)

    y = Variable((3, 3), :Semidefinite)
    p = minimize(trace(y), y[2,1]<=4, y[2,2]>=3)
    # @fact vexity(p) => ConvexVexity()
    solve!(p)
    @fact p.optval => roughly(3, TOL)

    x = Variable(Positive())
    y = Semidefinite(3)
    p = minimize(y[1, 2], y[2, 1] == 1)
    # @fact vexity(p) => ConvexVexity()
    solve!(p)
    @fact p.optval => roughly(1, TOL)
  end

  context("sdp constraints") do
    x = Variable(Positive())
    y = Variable((3, 3))
    p = minimize(x + y[1, 1], isposdef(y), x >= 1, y[2, 1] == 1)
    # @fact vexity(p) => ConvexVexity()
    solve!(p)
    @fact p.optval => roughly(1, TOL)
  end

  context("nuclear norm atom") do
    y = Semidefinite(3)
    p = minimize(nuclearnorm(y), y[2,1]<=4, y[2,2]>=3, y[3,3]<=2)
    @fact vexity(p) => ConvexVexity()
    solve!(p)
    @fact p.optval => roughly(3, TOL)
    @fact evaluate(nuclearnorm(y)) => roughly(3, TOL)
  end

  context("operator norm atom") do
    y = Variable((3,3))
    p = minimize(operatornorm(y), y[2,1]<=4, y[2,2]>=3, sum(y)>=12)
    @fact vexity(p) => ConvexVexity()
    solve!(p)
    @fact p.optval => roughly(4, TOL)
    @fact evaluate(operatornorm(y)) => roughly(4, TOL)
  end

  context("sigma max atom") do
    y = Variable((3,3))
    p = minimize(sigmamax(y), y[2,1]<=4, y[2,2]>=3, sum(y)>=12)
    @fact vexity(p) => ConvexVexity()
    solve!(p)
    @fact p.optval => roughly(4, TOL)
    @fact evaluate(sigmamax(y)) => roughly(4, TOL)
  end

  context("lambda max atom") do
    y = Semidefinite(3)
    p = minimize(lambdamax(y), y[1,1]>=4)
    @fact vexity(p) => ConvexVexity()
    solve!(p)
    @fact p.optval => roughly(4, TOL)
    @fact evaluate(lambdamax(y)) => roughly(4, TOL)
  end

  context("lambda min atom") do
    y = Semidefinite(3)
    p = maximize(lambdamin(y), trace(y)<=6)
    @fact vexity(p) => ConvexVexity()
    solve!(p)
    @fact p.optval => roughly(2, TOL)
    @fact evaluate(lambdamin(y)) => roughly(2, TOL)
  end

  context("matrix frac atom") do
    x = [1, 2, 3]
    P = Variable(3, 3)
    p = minimize(matrixfrac(x, P), P <= 2*eye(3), P >= 0.5 * eye(3))
    @fact vexity(p) => ConvexVexity()
    solve!(p)
    @fact p.optval => roughly(7, TOL)
    @fact evaluate(matrixfrac(x, P))[1] => roughly(7, TOL)
  end

  context("matrix frac atom both arguments variable") do
    x = Variable(3)
    P = Variable(3, 3)
    p = minimize(matrixfrac(x, P), lambdamax(P) <= 2, x[1] >= 1)
    @fact vexity(p) => ConvexVexity()
    solve!(p)
    @fact p.optval => roughly(.5, TOL)
    @fact evaluate(matrixfrac(x, P))[1] => roughly(.5, TOL)
  end

  context("sum largest eigs") do
    x = Semidefinite(3)
    p = minimize(sumlargesteigs(x, 2), x >= 1)
    solve!(p)
    @fact p.optval => roughly(3, TOL)
    @fact evaluate(x) => roughly(ones(3, 3), TOL)

    x = Semidefinite(3)
    p = minimize(sumlargesteigs(x, 2), [x[i,:] >= i for i=1:3]...)
    solve!(p)
    @fact p.optval => roughly(8.4853, TOL)

    x1 = Semidefinite(3)
    p1 = minimize(lambdamax(x1), x1[1,1]>=4)
    solve!(p1)

    x2 = Semidefinite(3)
    p2 = minimize(sumlargesteigs(x2, 1), x2[1,1]>=4)
    solve!(p2)

    @fact p1.optval => roughly(p2.optval, TOL)

    x1 = Semidefinite(3)
    p1 = minimize(lambdamax(x1), [x1[i,:] >= i for i=1:3]...)
    solve!(p1)

    x2 = Semidefinite(3)
    p2 = minimize(sumlargesteigs(x2, 1), [x2[i,:] >= i for i=1:3]...)
    solve!(p2)

    @fact p1.optval => roughly(p2.optval, TOL)
  end

  context("sgeomean_eig atom") do
    n = 5
    X = Variable(n,n)
    p = maximize(sgeomean_eig(X), eye(n) - X in :sdp); solve!(p)
    @fact vecnorm(X.value - eye(n)) => roughly(0, TOL*n^2)
  end

  context("equivalence of sgeomean_eig and logdet atoms") do
    #example is d-optimal design from CVX
    #web.cvxr.com/cvx/examples/cvxbook/Ch07_statistical_estim/html/expdesign.html
    M = 10
    angles1 = linspace(3*pi/4, pi, M)
    angles2 = linspace(0, -pi/2, M)
    V = [3.*cos(angles1)' 1.5.*cos(angles2)';
       3.*sin(angles1)' 1.5.*sin(angles2)']
    P = 2*M

    lam1 = Variable(P)
    p1 = maximize(sgeomean_eig(V * diagm(lam1) * V'), sum(lam1) == 1, lam1 >= 0)
    solve!(p1)

    @fact evaluate(sgeomean_eig(V * diagm(lam1) * V')) => roughly(p1.optval, TOL)
    @fact p1.optval => roughly(3.182, TOL)

    lam2 = Variable(P)
    p2 = maximize(logdet(V * diagm(lam2) * V'), sum(lam2) == 1, lam2 >= 0)
    solve!(p2)

    @fact evaluate(logdet(V * diagm(lam2) * V')) => roughly(p2.optval, TOL)
    @fact p2.optval => roughly(2.315, TOL)

    @fact lam1.value => roughly(lam2.value, TOL)
    @fact lam1.value[1] => roughly(0.5, TOL)
    @fact lam1.value[10] => roughly(0.5, TOL)
  end
end
