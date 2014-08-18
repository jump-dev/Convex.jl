include("../CVX.jl")
using CVX_refactor
x = Variable(Positive())
y = Variable((3,3),Semidefinite())
p = minimize(x+y[1,1])
solve!(p)