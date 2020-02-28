# # N queens

using Convex, GLPK, LinearAlgebra, SparseArrays, Test
aux(str) = joinpath(@__DIR__, "aux_files", str) # path to auxiliary files
include(aux("antidiag.jl"))

n = 8
# We encode the locations of the queens with a matrix of binary random variables.
x = Variable((n, n), :Bin)

# Now we impose the constraints: at most one queen on any anti-diagonal, at most one queen on any diagonal, and we must have exactly one queen per row and per column.
## At most one queen on any anti-diagonal
constr = Constraint[sum(antidiag(x, k)) <= 1 for k = -n+2:n-2]
## At most one queen on any diagonal
constr += Constraint[sum(diag(x, k)) <= 1 for k = -n+2:n-2]
## Exactly one queen per row and one queen per column
constr += Constraint[sum(x, dims=1) == 1, sum(x, dims=2) == 1]
p = satisfy(constr)
solve!(p, GLPK.Optimizer)

# Let us test the results:
for k = -n+2:n-2
	@test evaluate(sum(antidiag(x, k))) <= 1
	@test evaluate(sum(diag(x, k))) <= 1
end
@test all(evaluate(sum(x, dims=1)) .≈ 1)
@test all(evaluate(sum(x, dims=2)) .≈ 1)
