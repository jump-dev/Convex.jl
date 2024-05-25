# Convex.jl

[![CI](https://github.com/jump-dev/Convex.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/jump-dev/Convex.jl/actions/workflows/ci.yml)
[![Coverage](https://codecov.io/gh/jump-dev/Convex.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/jump-dev/Convex.jl)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://jump.dev/Convex.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://jump.dev/Convex.jl/dev)

Convex.jl is a [Julia](http://julialang.org) package for
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

## Getting help

For usage questions, please contact us via [Discourse](https://discourse.julialang.org/c/domain/opt).

If you have a reproducible example of a bug, please open a [GitHub issue](https://github.com/jump-dev/Convex.jl/issues/new).

## Installation

Install Convex using the Julia package manager:

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

# The problem is to minimize ||Ax - b||^2 subject to 0 <= x <= 1
# This can be done by: minimize(objective, constraints)
problem = minimize(sumsquares(A * x - b), [x >= 0, x <= 1])

# Solve the problem by calling solve!
solve!(problem, SCS.Optimizer)

# Check the status of the problem
problem.status

# Get the optimal value
problem.optval
```

## Using with JuMP

Convex.jl contains an experimental JuMP solver. This solver reformulates a
nonlinear JuMP model into a conic program using DCP. Note that it currently
supports only a limited subset of scalar nonlinear programs, such as those
involving `log` and `exp`.

```julia
julia> using JuMP, Convex, Clarabel

julia> model = Model(() -> Convex.Optimizer(Clarabel.Optimizer));

julia> set_silent(model)

julia> @variable(model, x >= 1);

julia> @variable(model, t);

julia> @constraint(model, t >= exp(x))
t - exp(x) ≥ 0

julia> @objective(model, Min, t);

julia> optimize!(model)

julia> value(x), value(t)
(0.9999999919393833, 2.7182818073461403)
```

## More Examples

A number of examples can be found [here](https://jump.dev/Convex.jl/stable/).
The [basic usage notebook](https://jump.dev/Convex.jl/stable/examples/general_examples/basic_usage/)
gives a simple tutorial on problems that can be solved using Convex.jl.

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
