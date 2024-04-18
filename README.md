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
termination_status(problem)

# Get the optimal value
objective_value(problem)
```

## Using with JuMP

The `master` branch of this package (not yet released) contains an experimental
JuMP solver. This solver reformulates a nonlinear JuMP model into a conic
program using DCP. Note that it currently supports only a limited subset of
scalar nonlinear programs, such as those involving `log` and `exp`.

```julia
julia> model = Model(() -> Convex.Optimizer(Clarabel.Optimizer))
A JuMP Model
Feasibility problem with:
Variables: 0
Model mode: AUTOMATIC
CachingOptimizer state: EMPTY_OPTIMIZER
Solver name: Convex with Clarabel

julia> @variable(model, x >= 1);

julia> @variable(model, t);

julia> @constraint(model, t >= exp(x))
t - exp(x) ≥ 0

julia> @objective(model, Min, t);

julia> optimize!(model)
-------------------------------------------------------------
           Clarabel.jl v0.5.1  -  Clever Acronym
                   (c) Paul Goulart
                University of Oxford, 2022
-------------------------------------------------------------

problem:
  variables     = 3
  constraints   = 5
  nnz(P)        = 0
  nnz(A)        = 5
  cones (total) = 2
    : Nonnegative = 1,  numel = 2
    : Exponential = 1,  numel = 3

settings:
  linear algebra: direct / qdldl, precision: Float64
  max iter = 200, time limit = Inf,  max step = 0.990
  tol_feas = 1.0e-08, tol_gap_abs = 1.0e-08, tol_gap_rel = 1.0e-08,
  static reg : on, ϵ1 = 1.0e-08, ϵ2 = 4.9e-32
  dynamic reg: on, ϵ = 1.0e-13, δ = 2.0e-07
  iter refine: on, reltol = 1.0e-13, abstol = 1.0e-12,
               max iter = 10, stop ratio = 5.0
  equilibrate: on, min_scale = 1.0e-04, max_scale = 1.0e+04
               max iter = 10

iter    pcost        dcost       gap       pres      dres      k/t        μ       step
---------------------------------------------------------------------------------------------
  0   0.0000e+00   4.4359e-01  4.44e-01  8.68e-01  8.16e-02  1.00e+00  1.00e+00   ------
  1   2.2037e+00   2.6563e+00  2.05e-01  7.34e-02  6.03e-03  5.44e-01  1.01e-01  9.33e-01
  2   2.5276e+00   2.6331e+00  4.17e-02  1.43e-02  1.26e-03  1.27e-01  2.26e-02  7.84e-01
  3   2.6758e+00   2.7129e+00  1.39e-02  4.09e-03  3.42e-04  4.35e-02  6.00e-03  7.84e-01
  4   2.7167e+00   2.7178e+00  3.90e-04  1.18e-04  9.82e-06  1.25e-03  1.72e-04  9.80e-01
  5   2.7182e+00   2.7183e+00  9.60e-06  3.39e-06  2.82e-07  3.15e-05  4.95e-06  9.80e-01
  6   2.7183e+00   2.7183e+00  1.92e-07  6.74e-08  5.62e-09  6.29e-07  9.84e-08  9.80e-01
  7   2.7183e+00   2.7183e+00  4.70e-09  1.94e-09  1.61e-10  1.59e-08  2.83e-09  9.80e-01
---------------------------------------------------------------------------------------------
Terminated with status = solved
solve time =  941μs

julia> value(x), value(t)
(0.9999999919393833, 2.7182818073461403)
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
