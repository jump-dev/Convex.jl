# Dualization

# Sometimes it can be much faster to solve the dual problem than the primal problem.
# Some solvers will automatically dualize the problem when heuristics deem it beneficial,
# and alternative DCP modeling languages like CVX will also automatically dualize the
# problem in some cases.
# Convex.jl does not automatically dualize any problem, but it is easy to manually do so
# with Dualization.jl.
# Here, we will solve a simple semidefinite program (from [issue #492](https://github.com/jump-dev/Convex.jl/issues/492)) to show how easy it is to dualize the problem,
# and that it can potentially provide speed ups.

# First we load our packages:
using LinearAlgebra
using Convex
using SCS
using Random
using Dualization

# Then we setup some test data.
Random.seed!(2022)
p = 50
Σ = Symmetric(randn(p, p))
Σ = Σ * Σ'

# Now we formulate and solve our primal problem:
d = Variable(p)
problem = maximize(sum(d), 0 ≤ d, d ≤ 1, Σ ⪰ Diagonal(d))
@elapsed solve!(problem, SCS.Optimizer; silent_solver = true)

# To solve the dual problem instead, we simply call `dual_optimizer` on our
# optimizer function:
@elapsed solve!(problem, dual_optimizer(SCS.Optimizer); silent_solver = true)
