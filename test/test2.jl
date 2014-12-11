using Base.Test
using Convex

TOL = 1e-3

# Addition, equality constraints, positive variables
x = Variable()
y = Variable(Positive())
p = minimize(x+y,x==1)
solve!(p)
@test_approx_eq_eps p.optval 1 TOL

# geq inequality constraints
x = Variable()
p = minimize(x, x>=1)
solve!(p)
@test_approx_eq_eps p.optval 1 TOL

# maximize, leq inequality constraints
x = Variable()
p = maximize(x, x<=2)
solve!(p)
@test_approx_eq_eps p.optval 2 TOL

# nonbinding geq inequality constraints
x = Variable()
p = maximize(x, x<=2, x>=1)
solve!(p)
@test_approx_eq_eps p.optval 2 TOL

# scalar multiplication
x = Variable()
p = minimize(2*x, x>=1)
solve!(p)
@test_approx_eq_eps p.optval 2 TOL

# indexing
x = Variable(2)
p = minimize(x[1], x>=1)
solve!(p)
@test_approx_eq_eps p.optval 1 TOL

# sums
x = Variable(2,2)
p = minimize(sum(x), x>=1)
solve!(p)
@test_approx_eq_eps p.optval 4 TOL

# sums
x = Variable(2,2)
p = minimize(sum(x) - 2*x[1,1], x>=1, x[1,1]<=2)
solve!(p)
@test_approx_eq_eps p.optval 1 TOL

# abs
x = Variable()
p = minimize(abs(x), x<=-1)
solve!(p)
@test_approx_eq_eps p.optval 1 TOL

x = Variable(2,2)
p = minimize(sum(abs(x)), x[2,2]>=1, x[1,1]>=1, x>=0)
solve!(p)
@test_approx_eq_eps p.optval 2 TOL

# Diag
x = Variable(2,2)
p = minimize(sum(diag(x,1)), x>=1)
solve!(p)
@test_approx_eq_eps p.optval 1 TOL

# .*
x = Variable(3, Positive())
p = maximize(sum(x.*[1,2,3]), x<=1)
solve!(p)
@test_approx_eq_eps p.optval 6 TOL

# ./
x = Variable(3, Positive())
p = maximize(sum(x./[1,2,3]), x<=1)
solve!(p)
@test_approx_eq_eps p.optval 11/6 TOL