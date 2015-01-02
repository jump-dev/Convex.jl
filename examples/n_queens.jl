using Convex, GLPKMathProgInterface
n = 8
x = Variable((n, n), :Bin)
# We can have at most one queen on any anti-diagonal
constr = Constraint[sum(antidiag(x, k)) <= 1 for k = -n+2:n-2]
# We can have at most one queen on any diagonal
constr += Constraint[sum(diag(x, k)) <= 1 for k = -n+2:n-2]
# We must have one queen per row and one queen per column
constr += Constraint[sum(x, 1) == 1, sum(x, 2) == 1]
p = satisfy(constr)
solve!(p, GLPKSolverMIP())

# Testing
for k = -n+2:n-2
	@assert evaluate(sum(antidiag(x, k))) <= 1
	@assert evaluate(sum(diag(x, k))) <= 1
end
@assert all(evaluate(sum(x, 1)) .== 1)
@assert all(evaluate(sum(x, 2)) .== 1)
