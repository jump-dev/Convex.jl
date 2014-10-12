using Base.Test
using Convex

TOL = 1e-2

# exp
y = Variable()
p = minimize(exp(y), y>=0)
println(conic_problem(p))
solve!(p, SCS.SCSMathProgModel())
@test_approx_eq_eps p.optval 1 TOL

# log
y = Variable()
p = maximize(log(y), y<=1)
println(conic_problem(p))
solve!(p, SCS.SCSMathProgModel())
@test_approx_eq_eps p.optval 0 TOL