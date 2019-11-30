# # Section Allocation
# Suppose you have $n$ students in a class who need to be assigned to $m$
# discussion sections. Each student needs to be assigned to exactly one section.
# Each discussion section should have between 6 and 10 students. Suppose an
# $n \times m$ preference matrix $P$ is given, where $P_{ij}$ gives student
# $i$'s ranking for section $j$ (1 would mean it is the student's top choice,
# 10,000 or a large number would mean the student can not attend that section).
#
# The goal will be to get an allocation matrix $X$, where $X_{ij} = 1$ if
# student $i$ is assigned to section $j$ and $0$ otherwise. 

using Convex, GLPKMathProgInterface
aux(str) = joinpath(@__DIR__, "aux", str) # path to auxiliary files

# Load our preference matrix, `P`
include(aux("data.jl"))

X = Variable(size(P), :Bin)

# We want every student to be assigned to exactly one section. So, every row
# must have exactly one non-zero entry. In other words, the sum of all the
# columns for every row is 1. We also want each section to have between 6 and 10
# students, so the sum of all the rows for every column should be between these.
constraints = [sum(X, dims=2) == 1, sum(X, dims=1) <= 10, sum(X, dims=1) >= 6]

# Our objective is simple `sum(X .* P)`, which can be more efficiently
# represented as `vec(X)' * vec(P)`. Since each entry of `X` is either 0 or 1,
# this is basically summing up the rankings of students that were assigned to them.
# If all students got their first choice, this value will be the number of
# students since the ranking of the first choice is 1.
p = minimize(vec(X)' * vec(P), constraints)

solve!(p, GLPKSolverMIP())
p.optval
