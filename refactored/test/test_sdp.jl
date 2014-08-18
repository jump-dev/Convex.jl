include("../CVX.jl")
using CVX_refactor

# SDP variables
x = Variable(Positive())
y = Variable((3,3), Semidefinite())
p = minimize(x+y[1,1])
solve!(p)


# SDP constraints
x = Variable(Positive())
y = Variable((3,3))
p = minimize(x + y[1,1], isposdef(y), x>=1,y[2,3]<=4)
solve!(p)

# nuclear norm
y = Variable((3,3), Semidefinite())
p = minimize(nuclear_norm(y), y[2,3]<=4, y[2,2]>=3)
solve!(p)

# trace
y = Variable((3,3), Semidefinite())
p = minimize(trace(y), y[2,3]<=4, y[2,2]>=3)
solve!(p)

# operator norm
y = Variable((3,3))
p = minimize(operator_norm(y), y[2,3]<=4, y[2,2]>=3, sum(y)>=12)
solve!(p)