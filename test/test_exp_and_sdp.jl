using Convex
using FactCheck

TOL = 1e-2

facts("SDP and Exp Atoms") do

  context("log det atom") do
    x = Variable(2, 2)
    p = maximize(logdet(x), [x[1, 1] == 1, x[2, 2] == 1])
    @fact vexity(p) --> ConvexVexity()
    solve!(p)
    @fact p.optval --> roughly(0, TOL)
    @fact evaluate(logdet(x)) --> roughly(0, TOL)
  end

end
