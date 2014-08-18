using Base.Test
include("../CVX.jl")
using CVX_refactor

TOL = 1e-3

# Addition, equality constraints, positive variables
x = Variable()
y = Variable(Positive())
p = minimize(x+y,x==1)
solve!(p)
@test_approx_eq_eps p.optval 1 TOL

# inequality constraints
x = Variable()
p = minimize(x, x>=1)
solve!(p)
@test_approx_eq_eps p.optval 1 TOL

# nonbinding inequality constraints
x = Variable()
p = maximize(x, x>=1, x<=2)
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
println(dual_conic_problem(p))
solve!(p)
@test_approx_eq_eps p.optval 4 TOL

# sums
x = Variable(2,2)
p = minimize(sum(x) - 2*x[1,1], x>=1, x[1,1]<=2)
println(dual_conic_problem(p))
solve!(p)
@test_approx_eq_eps p.optval 1 TOL

# Diag
x = Variable(2,2)
p = minimize(sum(diag(x,1)), x>=1)
solve!(p)
@test_approx_eq_eps x.value eye(2) TOL
