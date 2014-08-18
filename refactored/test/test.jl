include("../CVX.jl")
using CVX_refactor

# Addition, equality constraints, positive variables
x = Variable()
y = Variable(Positive())
p = minimize(x+y,x==1)
solve!(p)

# sums, inequality constraints, scalar mult, indexing
x = Variable(2,2)
p = minimize(sum(x) - 2*x[1,1], x>=1, x[1,1]<=2)
solve!(p)

# Diag
x = Variable(2,2)
p = minimize(sum(diag(x,1)), x>=1)

solve!(p)
