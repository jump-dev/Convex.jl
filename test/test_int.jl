using Convex
using FactCheck

TOL = 1e-2

LPsolver() = get_default_solver()

if !can_solve_mip(get_default_solver())
	using GLPKMathProgInterface
	MIPsolver() = GLPKSolverMIP()
else
	MIPsolver() = get_default_solver()
end

facts("Mixed Integer Programs") do

  context("lp fallback interface") do
    x = Variable()
    p = minimize(x, x>=4.3)
    @fact vexity(p) --> AffineVexity()
    solve!(p, LPsolver())
    @fact p.optval --> roughly(4.3, TOL)

    x = Variable(2)
    p = minimize(norm(x,1), x[1]>=4.3)
    @fact vexity(p) --> ConvexVexity()
    solve!(p, LPsolver())
    @fact p.optval --> roughly(4.3, TOL)
  end

  context("integer variables") do
    x = Variable(:Int)
    p = minimize(x, x>=4.3)
    @fact vexity(p) --> AffineVexity()
    solve!(p, MIPsolver())
    @fact p.optval --> roughly(5, TOL)

    x = Variable(2, :Int)
    p = minimize(sum(x), x>=4.3)
    @fact vexity(p) --> AffineVexity()
    solve!(p, MIPsolver())
    @fact p.optval --> roughly(10, TOL)

    x = Variable(:Int)
    y = Variable()
    p = minimize(sum(x + y), x>=4.3, y>=7)
    @fact vexity(p) --> AffineVexity()
    solve!(p, MIPsolver())
    @fact p.optval --> roughly(12, TOL)

    x = Variable(2, :Int)
    p = minimize(norm(x, 1), x[1]>=4.3)
    @fact vexity(p) --> ConvexVexity()
    solve!(p, MIPsolver())
    @fact p.optval --> roughly(5, TOL)

    x = Variable(2, :Int)
    p = minimize(sum(x), x[1]>=4.3, x>=0)
    @fact vexity(p) --> AffineVexity()
    solve!(p, MIPsolver())
    @fact p.optval --> roughly(5, TOL)

    x = Variable(2, :Int)
    p = minimize(sum(x), x>=.5)
    @fact vexity(p) --> AffineVexity()
    solve!(p, MIPsolver())
    @fact p.optval --> roughly(2, TOL)
  end

  context("binary variables") do
    x = Variable(2, :Bin)
    p = minimize(sum(x), x>=.5)
    @fact vexity(p) --> AffineVexity()
    solve!(p, MIPsolver())
    @fact p.optval --> roughly(2, TOL)

    x = Variable(2, :Bin)
    p = minimize(sum(x), x[1]>=.5, x>=0)
    @fact vexity(p) --> AffineVexity()
    solve!(p, MIPsolver())
    @fact p.optval --> roughly(1, TOL)
  end

end
