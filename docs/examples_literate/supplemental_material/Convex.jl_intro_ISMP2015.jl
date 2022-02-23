# # Convex Optimization in Julia
#
# ## Madeleine Udell | ISMP 2015

#-

# ## Convex.jl team
#
# * [Convex.jl](https://github.com/cvxgrp/Convex.jl): Madeleine Udell, Karanveer Mohan, David Zeng, Jenny Hong

#-

# ## Collaborators/Inspiration:
#
# * [CVX](http://www.cvxr.com): Michael Grant, Stephen Boyd
# * [CVXPY](https://github.com/cvxgrp/cvxpy): Steven Diamond, Eric Chu, Stephen Boyd
# * [JuliaOpt](https://github.com/JuliaOpt): Miles Lubin, Iain Dunning, Joey Huchette

## initial package installation

#-

## Make the Convex.jl module available
using Convex, SparseArrays, LinearAlgebra
using SCS # first order splitting conic solver [O'Donoghue et al., 2014]

## Generate random problem data
m = 50;  n = 100
A = randn(m, n)
x♮ = sprand(n, 1, .5) # true (sparse nonnegative) parameter vector
noise = .1*randn(m)    # gaussian noise
b = A*x♮ + noise      # noisy linear observations

## Create a (column vector) variable of size n.
x = Variable(n)

## nonnegative elastic net with regularization
λ = 1
μ = 1
problem = minimize(square(norm(A * x - b)) + λ*square(norm(x)) + μ*norm(x, 1),
                   x >= 0)

## Solve the problem by calling solve!
solve!(problem, MOI.OptimizerWithAttributes(SCS.Optimizer, "verbose" => 0))

println("problem status is ", problem.status) # :Optimal, :Infeasible, :Unbounded etc.
println("optimal value is ", problem.optval)

#-

using Interact, Plots
## Interact.WebIO.install_jupyter_nbextension() # might be helpful if you see `WebIO` warnings in Jupyter
@manipulate throttle=.1 for λ=0:.1:5, μ=0:.1:5
    global A
    problem = minimize(square(norm(A * x - b)) + λ*square(norm(x)) + μ*norm(x, 1),
                   x >= 0)
    solve!(problem, MOI.OptimizerWithAttributes(SCS.Optimizer, "verbose" => 0))
    histogram(evaluate(x), xlims=(0,3.5), label="x")
end

# # Quick convex prototyping

#-

# ## Variables

## Scalar variable
x = Variable()

#-

## (Column) vector variable
y = Variable(4)

#-

## Matrix variable
Z = Variable(4, 4)

# # Expressions

#-

# Convex.jl allows you to use a [wide variety of functions](http://convexjl.readthedocs.org/en/latest/operations.html) on variables and on expressions to form new expressions.

x + 2x

#-

e = y[1] + logdet(Z) + sqrt(x) + minimum(y)

# ### Examine the expression tree

e.children[2]

# # Constraints
#
# A constraint is convex if convex combinations of feasible points are also feasible. Equivalently, feasible sets are convex sets.
#
# In other words, convex constraints are of the form
#
# * `convexExpr <= 0`
# * `concaveExpr >= 0`
# * `affineExpr == 0`

x <= 0

#-

square(x) <= sum(y)

#-

M = Z
for i = 1:length(y)
    global M += rand(size(Z)...)*y[i]
end
M ⪰ 0

# # Problems

x = Variable()
y = Variable(4)
objective = 2*x + 1 - sqrt(sum(y))
constraint = x >= maximum(y)
p = minimize(objective, constraint)

#-

## solve the problem
solve!(p, MOI.OptimizerWithAttributes(SCS.Optimizer, "verbose" => 0))
p.status

#-

evaluate(x)

#-

## can evaluate expressions directly
evaluate(objective)

# ## Pass to solver
#
# call a `MathProgBase` solver suited for your problem class
#
# * see the [list of Convex.jl operations](http://convexjl.readthedocs.org/en/latest/operations.html) to find which cones you're using
# * see the [list of solvers](http://www.juliaopt.org/) for an up-to-date list of solvers and which cones they support

#-

# to solve problem using a different solver, just import the solver package and pass the solver to the `solve!` method: eg
#
#     using Mosek
#     solve!(p, Mosek.Optimizer)

#-

# ## Warmstart

## Generate random problem data
m = 50;  n = 100
A = randn(m, n)
x♮ = sprand(n, 1, .5) # true (sparse nonnegative) parameter vector
noise = .1*randn(m)    # gaussian noise
b = A*x♮ + noise      # noisy linear observations

## Create a (column vector) variable of size n.
x = Variable(n)

## nonnegative elastic net with regularization
λ = 1
μ = 1
problem = minimize(square(norm(A * x - b)) + λ*square(norm(x)) + μ*norm(x, 1),
                   x >= 0)
@time solve!(problem, MOI.OptimizerWithAttributes(SCS.Optimizer, "verbose" => 0))
λ = 1.5
@time solve!(problem, MOI.OptimizerWithAttributes(SCS.Optimizer, "verbose" => 0), warmstart = true)

# # DCP examples

## affine
x = Variable(4)
y = Variable(2)
sum(x) + y[2]

#-

2*maximum(x) + 4*sum(y) - sqrt(y[1] + x[1]) - 7 * minimum(x[2:4])

#-

## not dcp compliant
log(x) + square(x)

#-

## $f$ is convex increasing and $g$ is convex
square(pos(x))

#-

## $f$ is convex decreasing and $g$ is concave
invpos(sqrt(x))

#-

## $f$ is concave increasing and $g$ is concave
sqrt(sqrt(x))
