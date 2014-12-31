using Base.Test
using Convex

TOL = 1e-2

# SDP variables
y = Variable((2,2), :Semidefinite)
p = minimize(y[1,1])
solve!(p)
@test_approx_eq_eps p.optval 0 TOL

# SDP variables twice
y = Variable((3,3), :Semidefinite)
p = minimize(y[1,1], y[2,2]==1)
solve!(p)
@test_approx_eq_eps p.optval 0 TOL

# Solution is obtained as y[2,2] -> infinity
y = Variable((2, 2), :Semidefinite)
p = minimize(y[1, 1], y[1, 2] == 1)
solve!(p)
@test_approx_eq_eps p.optval 0 TOL

# SDP variables
y = Semidefinite(3)
p = minimize(sum(diag(y)), y[1, 1] == 1)
solve!(p)
@test_approx_eq_eps p.optval 1 TOL

# SDP constraints
x = Variable(Positive())
y = Variable((3, 3))
p = minimize(x + y[1, 1], isposdef(y), x >= 1, y[2, 1] == 1)
solve!(p)
@test_approx_eq_eps p.optval 1 TOL

x = Variable(Positive())
y = Semidefinite(3)
p = minimize(y[1, 2], y[2, 1] == 1)
solve!(p)
@test_approx_eq_eps p.optval 1 TOL

# Not symmetric
x = Variable(Positive())
y = Semidefinite(3, is_symmetric=false)
p = minimize(y[1, 2], y[2, 1] == 1, y[1, 2] >= -1000)
solve!(p)
@test_approx_eq_eps p.optval -1000 TOL

# trace
y = Variable((3, 3), :Semidefinite)
p = minimize(trace(y), y[2,1]<=4, y[2,2]>=3)
solve!(p)
@test_approx_eq_eps p.optval 3 TOL

# nuclear norm
y = Semidefinite(3)
p = minimize(nuclear_norm(y), y[2,1]<=4, y[2,2]>=3, y[3,3]<=2)
solve!(p)
@test_approx_eq_eps p.optval 3 TOL

# operator norm
y = Variable((3,3))
p = minimize(operator_norm(y), y[2,1]<=4, y[2,2]>=3, sum(y)>=12)
solve!(p)
@test_approx_eq_eps p.optval 4 TOL

# sigma_max
y = Variable((3,3))
p = minimize(sigma_max(y), y[2,1]<=4, y[2,2]>=3, sum(y)>=12)
solve!(p)
@test_approx_eq_eps p.optval 4 TOL

# lambda_max
y = Semidefinite(3)
p = minimize(lambda_max(y), y[1,1]>=4)
solve!(p)
@test_approx_eq_eps p.optval 4 TOL

# lambda_min
y = Semidefinite(3)
p = maximize(lambda_min(y), trace(y)<=6)
solve!(p)
@test_approx_eq_eps p.optval 2 TOL