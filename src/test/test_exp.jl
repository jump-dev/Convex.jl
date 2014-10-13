using Base.Test
using Convex

TOL = 1e-2

# exp
y = Variable()
p = minimize(exp(y), y>=0)
solve!(p, SCS.SCSMathProgModel())
@test_approx_eq_eps p.optval 1 TOL

y = Variable()
p = minimize(exp(y), y>=1)
solve!(p, SCS.SCSMathProgModel())
@test_approx_eq_eps p.optval exp(1) TOL

# log
y = Variable()
p = maximize(log(y), y<=1)
solve!(p, SCS.SCSMathProgModel())
@test_approx_eq_eps p.optval 0 TOL

y = Variable()
p = maximize(log(y), y<=2)
solve!(p, SCS.SCSMathProgModel())
@test_approx_eq_eps p.optval log(2) TOL

y = Variable()
p = maximize(log(y), [y<=2, exp(y)<=10])
solve!(p, SCS.SCSMathProgModel())
@test_approx_eq_eps p.optval log(2) TOL

# multidimensional exp
y = Variable(5)
p = minimize(sum(exp(y)), y>=0)
solve!(p, SCS.SCSMathProgModel())
@test_approx_eq_eps p.optval 5 TOL