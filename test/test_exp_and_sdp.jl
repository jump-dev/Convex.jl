using Base.Test
using Convex
# using SCS
# set_default_solver(SCSSolver())

TOL = 1e-2

x = Variable(2, 2)
p = maximize(logdet(x), [x[1, 1] == 1, x[2, 2] == 1])
solve!(p)
@test_approx_eq_eps p.optval 0 TOL
