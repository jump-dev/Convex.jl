# Convex.jl

[![Build Status](https://github.com/jump-dev/Convex.jl/workflows/CI/badge.svg)](https://github.com/jump-dev/Convex.jl/actions?query=workflow%3ACI)
[![Coverage](https://codecov.io/gh/jump-dev/Convex.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/jump-dev/Convex.jl)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://jump.dev/Convex.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://jump.dev/Convex.jl/dev)

**Convex.jl** is a [Julia](http://julialang.org) package for
[Disciplined Convex Programming](http://dcp.stanford.edu/) (DCP).

Convex.jl can solve linear programs, mixed-integer linear programs, and
DCP-compliant convex programs using a variety of solvers, including
[Mosek](https://github.com/MOSEK/Mosek.jl),
[Gurobi](https://github.com/jump-dev/Gurobi.jl),
[ECOS](https://github.com/jump-dev/ECOS.jl),
[SCS](https://github.com/jump-dev/SCS.jl), and
[GLPK](https://github.com/jump-dev/GLPK.jl), through
[MathOptInterface](https://github.com/jump-dev/MathOptInterface.jl).

Convex.jl also supports optimization with complex variables and coefficients.

For usage questions, please contact us via [Discourse](https://discourse.julialang.org/c/domain/opt).

## Installation

```julia
import Pkg
Pkg.add("Convex")
```

## Quick Example

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
gives a simple tutorial on problems that can be solved using Convex.jl.

All examples can be downloaded as a zip file from [here](https://jump.dev/Convex.jl/stable/examples/notebooks.zip).

## Citing this package

If you use Convex.jl for published work, we encourage you to cite the software
using the following BibTeX citation:

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
