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
    p = minimize(nuclear_norm(y), y[2,1]<=4, y[2,2]>=3, y[3,3]<=2)
    @fact vexity(p) => ConvexVexity()
    solve!(p)
    @fact p.optval => roughly(3, TOL)
    @fact evaluate(nuclear_norm(y)) => roughly(3, TOL)
  end

  context("operator norm atom") do
    y = Variable((3,3))
    p = minimize(operator_norm(y), y[2,1]<=4, y[2,2]>=3, sum(y)>=12)
    @fact vexity(p) => ConvexVexity()
    solve!(p)
    @fact p.optval => roughly(4, TOL)
    @fact evaluate(operator_norm(y)) => roughly(4, TOL)
  end

  context("sigma max atom") do
    y = Variable((3,3))
    p = minimize(sigma_max(y), y[2,1]<=4, y[2,2]>=3, sum(y)>=12)
    @fact vexity(p) => ConvexVexity()
    solve!(p)
    @fact p.optval => roughly(4, TOL)
    @fact evaluate(sigma_max(y)) => roughly(4, TOL)
  end

  context("lambda max atom") do
    y = Semidefinite(3)
    p = minimize(lambda_max(y), y[1,1]>=4)
    @fact vexity(p) => ConvexVexity()
    solve!(p)
    @fact p.optval => roughly(4, TOL)
    @fact evaluate(lambda_max(y)) => roughly(4, TOL)
  end

  context("lambda min atom") do
    y = Semidefinite(3)
    p = maximize(lambda_min(y), trace(y)<=6)
    @fact vexity(p) => ConvexVexity()
    solve!(p)
    @fact p.optval => roughly(2, TOL)
    @fact evaluate(lambda_min(y)) => roughly(2, TOL)
  end

  context("matrix frac atom") do
    x = [1, 2, 3]
    P = Variable(3, 3)
    p = minimize(matrix_frac(x, P), P <= 2*eye(3), P >= 0.5 * eye(3))
    @fact vexity(p) => ConvexVexity()
    solve!(p)
    @fact p.optval => roughly(7, TOL)
    @fact evaluate(matrix_frac(x, P))[1] => roughly(7, TOL)
  end
end
