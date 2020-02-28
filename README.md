# Convex.jl

[![Build Status](https://travis-ci.org/JuliaOpt/Convex.jl.svg?branch=master)](https://travis-ci.org/JuliaOpt/Convex.jl)
[![Coverage Status](https://coveralls.io/repos/JuliaOpt/Convex.jl/badge.svg?branch=master)](https://coveralls.io/r/JuliaOpt/Convex.jl?branch=master)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://www.juliaopt.org/Convex.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://www.juliaopt.org/Convex.jl/dev)

**Convex.jl** is a [Julia](http://julialang.org) package for [Disciplined Convex Programming](http://dcp.stanford.edu/). Convex.jl can solve linear programs, mixed-integer linear programs, and DCP-compliant convex programs using a variety of solvers, including [Mosek](https://github.com/JuliaOpt/Mosek.jl), [Gurobi](https://github.com/JuliaOpt/Gurobi.jl), [ECOS](https://github.com/JuliaOpt/ECOS.jl), [SCS](https://github.com/JuliaOpt/SCS.jl), and  [GLPK](https://github.com/JuliaOpt/GLPK.jl), through [MathOptInterface](https://github.com/JuliaOpt/MathOptInterface.jl). It also supports optimization with complex variables and coefficients.

**Installation**: `julia> Pkg.add("Convex")`

- **Detailed documentation and examples** for Convex.jl ([stable](https://juliaopt.github.io/Convex.jl/stable) | [development version](https://juliaopt.github.io/Convex.jl/dev)).
- If you're running into **bugs or have feature requests**, please use the [Github Issue Tracker](https://github.com/JuliaOpt/Convex.jl/issues>).
- For usage questions, please contact us via [Discourse](https://discourse.julialang.org/c/domain/opt).

## Quick Example

Note: the following example uses the syntax of the master branch of Convex.jl. See the [stable docs](https://www.juliaopt.org/Convex.jl/stable) for examples for the current release.

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

A number of examples can be found [here](https://www.juliaopt.org/Convex.jl/stable/).
The [basic usage notebook](https://www.juliaopt.org/Convex.jl/stable/examples/general_examples/basic_usage/)
gives a simple tutorial on problems that can be solved using Convex.jl. All examples can be downloaded as
a zip file from [here](https://www.juliaopt.org/Convex.jl/stable/examples/notebooks.zip).

## Citing this package

If you use Convex.jl for published work, we encourage you to cite the software using the following BibTeX citation:
```
@article{convexjl,
 title = {Convex Optimization in {J}ulia},
 author ={Udell, Madeleine and Mohan, Karanveer and Zeng, David and Hong, Jenny and Diamond, Steven and Boyd, Stephen},
 year = {2014},
 journal = {SC14 Workshop on High Performance Technical Computing in Dynamic Languages},
 archivePrefix = "arXiv",
 eprint = {1410.4821},
 primaryClass = "math-oc",
}
```

> Convex.jl was previously called CVX.jl.
