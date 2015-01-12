using Convex
using FactCheck

TOL = 1e-2

facts("SDP Atoms") do
  context("sdp variables") do
    y = Variable((2,2), :Semidefinite)
    p = minimize(y[1,1])
    solve!(p)
    @fact p.optval => roughly(0, TOL)

    y = Variable((3,3), :Semidefinite)
    p = minimize(y[1,1], y[2,2]==1)
    solve!(p)
    @fact p.optval => roughly(0, TOL)

    # Solution is obtained as y[2,2] -> infinity
    # This test fails on Mosek. See
    # https://github.com/JuliaOpt/Mosek.jl/issues/29
    y = Variable((2, 2), :Semidefinite)
    p = minimize(y[1, 1], y[1, 2] == 1)
    solve!(p)
    @fact p.optval => roughly(0, TOL)

    y = Semidefinite(3)
    p = minimize(sum(diag(y)), y[1, 1] == 1)
    solve!(p)
    @fact p.optval => roughly(1, TOL)

    y = Variable((3, 3), :Semidefinite)
    p = minimize(trace(y), y[2,1]<=4, y[2,2]>=3)
    solve!(p)
    @fact p.optval => roughly(3, TOL)

    x = Variable(Positive())
    y = Semidefinite(3)
    p = minimize(y[1, 2], y[2, 1] == 1)
    solve!(p)
    @fact p.optval => roughly(1, TOL)
  end

  context("sdp constraints") do
    x = Variable(Positive())
    y = Variable((3, 3))
    p = minimize(x + y[1, 1], isposdef(y), x >= 1, y[2, 1] == 1)
    solve!(p)
    @fact p.optval => roughly(1, TOL)
  end

  context("nuclear norm atom") do
    y = Semidefinite(3)
    p = minimize(nuclear_norm(y), y[2,1]<=4, y[2,2]>=3, y[3,3]<=2)
    solve!(p)
    @fact p.optval => roughly(3, TOL)
  end

  context("operator norm atom") do
    y = Variable((3,3))
    p = minimize(operator_norm(y), y[2,1]<=4, y[2,2]>=3, sum(y)>=12)
    solve!(p)
    @fact p.optval => roughly(4, TOL)
  end

  context("sigma max atom") do
    y = Variable((3,3))
    p = minimize(sigma_max(y), y[2,1]<=4, y[2,2]>=3, sum(y)>=12)
    solve!(p)
    @fact p.optval => roughly(4, TOL)
  end

  context("lambda max atom") do
    y = Semidefinite(3)
    p = minimize(lambda_max(y), y[1,1]>=4)
    solve!(p)
    @fact p.optval => roughly(4, TOL)
  end

  context("lambda min atom") do
    y = Semidefinite(3)
    p = maximize(lambda_min(y), trace(y)<=6)
    solve!(p)
    @fact p.optval => roughly(2, TOL)
  end

end
