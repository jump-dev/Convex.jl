x# Please read: Convex.jl needs a new maintainer

Hi there! Convex.jl is an old Julia package that has had a number of maintainers
and [contributors over the years](https://github.com/jump-dev/Convex.jl/graphs/contributors).
Most typically, these have been procrastinating graduate students who work on
Convex.jl instead of finishing their thesis. However, we can't all be students
forever, and so is Convex.jl in need of a new maintainer. 

Therefore, **Convex.jl is in maintenance mode.**

There are still a number of community members who can review and merge pull
requests, but they will not be doing any active development.

If you would like to get involved and help revive Convex.jl:
 * Join the JuMP developer chatroom at [https://gitter.im/JuliaOpt/JuMP-dev](https://gitter.im/JuliaOpt/JuMP-dev)
   to say hi. If you have questions relating to Convex.jl, we may be able to help.
 * Start by making some small pull requests to fix any of the open issues.

# Convex.jl

[![Build Status](https://github.com/jump-dev/Convex.jl/workflows/CI/badge.svg)](https://github.com/jump-dev/Convex.jl/actions?query=workflow%3ACI)
[![Coverage](https://codecov.io/gh/jump-dev/Convex.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/jump-dev/Convex.jl)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://jump.dev/Convex.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://jump.dev/Convex.jl/dev)

**Convex.jl** is a [Julia](http://julialang.org) package for [Disciplined Convex Programming](http://dcp.stanford.edu/). Convex.jl can solve linear programs, mixed-integer linear programs, and DCP-compliant convex programs using a variety of solvers, including [Mosek](https://github.com/JuliaOpt/Mosek.jl), [Gurobi](https://github.com/jump-dev/Gurobi.jl), [ECOS](https://github.com/jump-dev/ECOS.jl), [SCS](https://github.com/jump-dev/SCS.jl), and  [GLPK](https://github.com/JuliaOpt/GLPK.jl), through [MathOptInterface](https://github.com/jump-dev/MathOptInterface.jl). It also supports optimization with complex variables and coefficients.

**Installation**: `julia> Pkg.add("Convex")`

- **Detailed documentation and examples** for Convex.jl ([stable](https://jump.dev/Convex.jl/stable) | [development version](https://jump.dev/Convex.jl/dev)).
- If you're running into **bugs or have feature requests**, please use the [Github Issue Tracker](https://github.com/jump-dev/Convex.jl/issues>).
- For usage questions, please contact us via [Discourse](https://discourse.julialang.org/c/domain/opt).

## Quick Example

To run this example, first install Convex and at least one solver, such as SCS:
```julia
using Pkg
Pkg.add("Convex")
Pkg.add("SCS")
```
Now let's solve a least-squares problem with inequality constraints.
```julia
# Let us first make the Convex.jl module available
using Convex, SCS

# Generate random problem data
m = 4;  n = 5
A = randn(m, n); b = randn(m, 1)

# Create a (column vector) variable of size n x 1.
x = Variable(n)

# The problem is to minimize ||Ax - b||^2 subject to x >= 0
# This can be done by: minimize(objective, constraints)
problem = minimize(sumsquares(A * x - b), [x >= 0])

# Solve the problem by calling solve!
solve!(problem, SCS.Optimizer)

# Check the status of the problem
problem.status # :Optimal, :Infeasible, :Unbounded etc.

# Get the optimal value
problem.optval
```

## More Examples

A number of examples can be found [here](https://jump.dev/Convex.jl/stable/).
The [basic usage notebook](https://jump.dev/Convex.jl/stable/examples/general_examples/basic_usage/)
gives a simple tutorial on problems that can be solved using Convex.jl. All examples can be downloaded as
a zip file from [here](https://jump.dev/Convex.jl/stable/examples/notebooks.zip).

## Citing this package

If you use Convex.jl for published work, we encourage you to cite the software using the following BibTeX citation:
```bibtex
@article{convexjl,
 title = {Convex Optimization in {J}ulia},
 author = {Udell, Madeleine and Mohan, Karanveer and Zeng, David and Hong, Jenny and Diamond, Steven and Boyd, Stephen},
 year = {2014},
 journal = {SC14 Workshop on High Performance Technical Computing in Dynamic Languages},
 archivePrefix = "arXiv",
 eprint = {1410.4821},
 primaryClass = "math-oc",
}
```

> Convex.jl was previously called CVX.jl.
