using Base.Test
using Convex
using GLPKMathProgInterface

TOL = 1e-2

MIPsolver() = GLPKSolverMIP()
LPsolver() = GLPKSolverLP()

# LP fallback interface
x = Variable()
p = minimize(x, x>=4.3)
solve!(p, LPsolver())
@test_approx_eq_eps p.optval 4.3 TOL

x = Variable(2)
p = minimize(norm(x,1), x[1]>=4.3)
solve!(p, LPsolver())
@test_approx_eq_eps p.optval 4.3 TOL

# integer variables
x = Variable(:Int)
p = minimize(x, x>=4.3)
solve!(p, MIPsolver())
@test_approx_eq_eps p.optval 5 TOL

x = Variable(2, :Int)
p = minimize(sum(x), x>=4.3)
solve!(p, MIPsolver())
@test_approx_eq_eps p.optval 10 TOL

x = Variable(:Int)
y = Variable()
p = minimize(sum(x + y), x>=4.3, y>=7)
solve!(p, MIPsolver())
@test_approx_eq_eps p.optval 12 TOL

x = Variable(2, :Int)
p = minimize(norm(x, 1), x[1]>=4.3)
solve!(p, MIPsolver())
@test_approx_eq_eps p.optval 5 TOL

x = Variable(2, :Int)
p = minimize(sum(x), x[1]>=4.3, x>=0)
solve!(p, MIPsolver())
@test_approx_eq_eps p.optval 5 TOL

x = Variable(2, :Int)
p = minimize(sum(x), x>=.5)
solve!(p, MIPsolver())
@test_approx_eq_eps p.optval 2 TOL

x = Variable(2, :Bin)
p = minimize(sum(x), x>=.5)
solve!(p, MIPsolver())
@test_approx_eq_eps p.optval 2 TOL

x = Variable(2, :Bin)
p = minimize(sum(x), x[1]>=.5, x>=0)
solve!(p, MIPsolver())
@test_approx_eq_eps p.optval 1 TOL
