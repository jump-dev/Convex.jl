# Convex.jl

[![Build Status](https://travis-ci.org/cvxgrp/Convex.jl.svg?branch=master)](https://travis-ci.org/cvxgrp/Convex.jl) [![Coverage Status](https://img.shields.io/coveralls/cvxgrp/Convex.jl.svg)](https://coveralls.io/r/cvxgrp/Convex.jl)

Convex.jl is a julia package for [Disciplined Convex Programming](http://dcp.stanford.edu/). Convex.jl can solve linear programs, mixed-integer linear programs, and dcp-compliant convex programs using a variety of solvers, including [Mosek](https://github.com/JuliaOpt/Mosek.jl), [Gurobi](https://github.com/JuliaOpt/gurobi.jl), [ECOS](https://github.com/JuliaOpt/ECOS.jl), [SCS](https://github.com/karanveerm/SCS.jl), [GLPK](https://github.com/JuliaOpt/GLPK.jl), through the [MathProgBase](http://mathprogbasejl.readthedocs.org/en/latest/) interface.

Note that Convex.jl was previously called CVX.jl. This package is under active development; we welcome bug reports and feature requests. For usage questions, please contact us via the [JuliaOpt mailing list](https://groups.google.com/forum/#!forum/julia-opt).

## Table of Contents
- [Introduction](#user-content-introduction)
- [Installation](#user-content-installation)
- [Basic Types](#user-content-basic-types)
  - [Variables](#user-content-variables)
  - [Constants](#user-content-constants)
  - [Expressions](#user-content-expressions)
  - [Constraints](#user-content-constraints)
  - [Objective](#user-content-objective)
  - [Problem](#user-content-problem)
- [Supported Operations](#user-content-supported-operations)
- [Examples](#user-content-examples)
- [Solvers](#user-content-solvers)
- [Credits](#user-content-credits)
- [Citing this package](#user-content-citing-this-package)

## Introduction
Convex.jl allows you to model and solve optimization problems by expressing them
in a simple, mathematical form. It is compatible with any solver in MathProgBase,
and can solve LPs, MIPs, SOCPs, SDPs, and exponential cone programs,
often, without the user needing to know what these are.
If your problem can be expressed while following the rules of
[Disciplined Convex Programming](http://dcp.stanford.edu/),
Convex.jl will transform your problem into a standard form that can be solved
using a solver of your choice.
For a detailed discussion of how Convex.jl works, see
[our paper](http://www.arxiv.org/abs/1410.4821).

Here's a quick example of code that solves a least-squares problem with inequality constraints:
```
using Convex

# Generate problem data
m = 4
n = 5
A = randn(m, n)
b = randn(m, 1)

# Create a (column vector) variable of size n x 1.
x = Variable(n)

# The problem is to minimize ||Ax + b||^2 subject to x >= 0
problem = minimize(sum_squares(A * x + b), [x >= 0])

# We can also add constraints after the problem is created in the following way:
problem.constraints += [x <= 1, 0.5 <= 2*x]

# Solve the problem by calling solve!
solve!(problem)

# Alternatively, specify a solver yourself
using ECOS
solver = ECOSSolver() # or SCSSolver() or GurobiSolver() or GLPKSolverMIP() etc.
solve!(problem, solver)

# Status (optimal, infeasible, etc.)
problem.status

# Optimal value
problem.optval

# Optimal value of x
x.value
```

## Installation
Convex.jl can be installed using the command
```
Pkg.add("Convex")
```
Only the solver ECOS is installed by default.
You can also add and use any solver in [JuliaOpt](https://github.com/JuliaOpt),
including GLPK, CPLEX, Gurobi, and MOSEK,
for LPs or MILPs. For more information on solvers see the [Solvers](#user-content-solvers) section below,
or refer to the documentation for that solver for installation instructions.

Currently [SCS](www.github.com/karanveerm/SCS.jl) is the only solver that can be used
to solve SDPs and exponential cone programs,
and SCS currently works only on OSX, so SDPs and exponential cone programs
are only supported on OSX for now. SCS can be installed using the following commands:
```
Pkg.clone("https://github.com/karanveerm/SCS.jl.git")
Pkg.build("SCS")
```
You might have to restart Julia.
Please file an issue in case you run into problems during installation. We'll be glad to help!

## Basic Types

The basic building block of Convex.jl is called an *expression*,
which can represent a variable, a constant, or a function of another expression.
We discuss each kind of expression in turn.

### Variables
The simplest kind of expression in Convex.jl is a variable.
Variables in Convex.jl are declared using the `Variable`
keyword, along with the dimensions of the variable.
```
# Scalar variable
x = Variable()

# Column-vector variable
x = Variable(5)

# Matrix variable
x = Variable(4, 6)
```

Variables may also be declared as having special properties, such as being

  * (entrywise) positive: `x = Variable(4, Positive())`,
  * (entrywise) negative: `x = Variable(4, Negative())`,
  * integral: `x = Variable(4, :Int)`,
  * binary: `x = Variable(4, :Bin)`, or
  * (for a matrix) being symmetric, with nonnegative eigenvalues (ie, positive semidefinite): `z = Semidefinite(4)`.

### Constants
Numbers, vectors, and matrices present in the Julia environment are wrapped
automatically into a `Constant` expression when used in a Convex.jl expression.

### Expressions
Expressions in Convex.jl are formed by applying any *atom* (mathematical function
defined in Convex.jl) to variables, constants, and other expressions.
For a list of these functions, see [Supported Operations](#user-content-supported-operations) below.
Atoms are applied to expressions using operator overloading. Hence, `2+2`
calls Julia's built-in addition operator, while `2+x` calls the Convex.jl
addition method and returns a Convex.jl expression.
Many of the useful language
features in Julia, such as arithmetic, array indexing, and matrix transpose are
overloaded in Convex.jl so they may be used with variables and expressions
just as they are used with native Julia types.

Expressions that are created must be DCP-compliant.
More information on DCP can be found [here](http://dcp.stanford.edu/).
```
x = Variable(5)
# The following are all expressions
y = sum(x)
z = 4 * x + y
z_1 = z[1]
```
Convex.jl allows the values of the expressions to be evaluated directly.
```
x = Variable()
y = Variable()
z = Variable()
expr = x + y + z
problem = minimize(expr, x >= 1, y >= x, 4 * z >= y)
solve!(problem)

# Once the problem is solved, we can call evaluate() on expr:
evaluate(expr)
```

### Constraints
*Constraints* in Convex.jl are declared using the standard comparison
operators `<=`, `>=`, and `==`.  They specify relations that
must hold between two expressions.  Convex.jl does not distinguish between strict
and non-strict inequality constraints.
```
x = Variable(5, 5)
# Equality constraint
constraint = x == 0
# Inequality constraint
constraint = x >= 1
```
Matrices can also be constrained to be positive semidefinite.
```
x = Variable(3, 3)
y = Variable(3, 1)
z = Variable()
# constrain [x y; y' z] to be positive semidefinite
constraint = isposdef([x y; y' z])
```

### Objective
The objective of the problem is a scalar expression to be maximized or minimized by using `maximize` or `minimize` respectively. Feasibility problems are also allowed by either giving a constant as the expression, or using `problem = satisfy(constraints)`.

### Problem
A *problem* in Convex.jl consists of a *sense* (minimize, maximize,
or satisfy), an objective (an expression to which the sense verb is to be
applied), and zero or more constraints which must be satisfied at the
solution.
Problems may be constructed as
```problem = minimize(objective, constraints)```
or
```problem = maximize(objective, constraints)```
or
```problem = satisfy(constraints)```

Constraints can be added at any time before the problem is solved.
```
# No constraints given
problem = minimize(objective)
# Add some constraint
problem.constraints += constraint
# Add many more constraints
problem.constraints += [constraint1, constraint2]
```
A problem can be solved by calling `solve!`:
```
solve!(problem)
```
After the problem is solved, `problem.status` records the status returned by the optimization solver,
and can be `:Optimal`, `:Infeasible`, `:Unbounded`, `:Indeterminate` or `:Error`.
If the status is `:Optimal`, `problem.optval` will record the optimum value of the problem.
The optimal value for each variable `x` participating in the problem
can be found in `x.value`.
The optimal value of an expression can be found by calling the `evaluate()` function
on the expression as follows: `evaluate(expr)`.
<!--The dual values are stored with the respective constraints and can be accessed as `problem.constraints[idx].dual_value`.-->

## Supported operations
Convex.jl currently supports the following operations.
These functions ("atoms") may be composed according to the
[DCP](dcp.stanford.edu) composition rules to form new convex, concave, or affine expressions.

### Affine atoms

These atoms are affine in their arguments.

 - addition, subtraction, multiplication, division: `+, -, /, *`
 - indexing into vectors and matrices: `x[1:4, 2:3]`
 - k-th diagonal of a matrix: `diag(x, k)`
 - transpose: `x'`
 - dot product: `x' * y` or `dot(x, y)`
 - reshape, vec: `reshape(x, 2, 3)` or `vec(x)`
 - minimum, maximum element of a vector or matrix: `maximum(x)`
 - horizontal and vertical stacking: `hcat(x, y); vcat(x, y)`

### Elementwise atoms

These atoms are applied elementwise to their arguments, returning an expression of the same size.

atom | description | vexity | slope | implicit constraint
-----|--------|--------|------- |-------
`min(x,y)` | $min(x,y)$ | convex | increasing | none
`max(x,y)` | $min(x,y)$ | convex | increasing | none
`pos(x)` | $max(x,0)$ | convex | increasing | none
`neg(x)` | $max(-x,0)$ | convex | decreasing | none
`inv_pos(x)` | $1/max(x,0)$ | convex | decreasing | $x>0$
`sqrt(x)` | $sqrt(x)$ | convex | decreasing | $x>0$
`square(x)`, `x^2` | $x^2$ | convex | increasing on $x \ge 0$, decreasing on $x \le 0$ | none
`abs(x)` | $abs(x)$ | convex | increasing on $x \ge 0$, decreasing on $x \le 0$ | none
`geo_mean(x, y)` | $\sqrt{xy}$ | concave | increasing | $x \ge 0$, $y \ge 0$
`exp(x)` | $\exp(x)$ | convex | increasing | none
`log(x)` | $\log(x)$ | concave | increasing | $x \gt 0$

### Vector and Matrix atoms

These atoms take vector or matrix arguments and return scalar expressions.

atom | description | vexity | slope | implicit constraint
-----|--------|--------|------- |-------
`norm(x, p)` | $(\sum x_i^p)^{1/p}$ | convex | increasing on $x \ge 0$, decreasing on $x \le 0$ | `p = 1, 2, Inf`
`vecnorm(x, p)` | $(\sum x_{ij}^p)^{1/p}$ | convex | increasing on $x \ge 0$, decreasing on $x \le 0$ | `p = 1, 2, Inf`
`quad_form(x,P)` | $x^T P x$ | convex in $X$, affine in $P$ | increasing on $x \ge 0$, decreasing on $x \le 0$, increasing in $P$ | either $x$ or $P$ must be constant
`quad_over_lin(x, y)` | $x^T x/y$ | convex | increasing on $x \ge 0$, decreasing on $x \le 0$, decreasing in $y$ | $y > 0$
`sum_squares(x)` | $\sum x_i^2$ | convex | increasing on $x \ge 0$, decreasing on $x \le 0$ | none
`nuclear_norm(x)` | sum of singular values of $x$ | convex | not monotonic | none
`operator_norm(x)` | maximum singular values of $x$ | convex | not monotonic | none
`logsumexp(x)` | $\log(\sum_i \exp(x_i))$ | convex | increasing | none


### Promotion
When an atom or constraint is applied to a scalar and a higher dimensional variable, the scalars are promoted. For example, we can do `max(x, 0)` gives an expression with the shape of `x` whose elements are the maximum of the corresponding element of `x` and `0`.

## Examples
A number of very simple examples can be found in the test/ directory.
More sophisticated examples, along with plots can be found in the examples/ directory.
Here are a few simple examples to start with:

* Dot Product
```
x = Variable(2)
A = 1.5 * eye(2)
p = minimize(dot([2.0; 2.0], x), A * x >= [1.1; 1.1])
solve!(p)
println(p.optval)
println(x.value)
```

* Matrix Variables
```
X = Variable(2, 2)
c = ones(2, 1)
p = minimize(c' * X * c, X >= ones(2, 2))
solve!(p)
println(X.value)
```

* Promotion of scalar variables and constants when present with matrices
```
N = 20
x = Variable(1) # x is a scalar variable
y = Variable(N, N) # y is a 20 x 20 matrix variavble
c = ones(N, 1)
# We can add scalar variables to matrix variables.
# The following line is equivalent to c' * (y + eye(N) * x) * c
# Similar overloading is done in other cases
objective = c' * (y + x) * c
p = minimize(objective, x >= 3, 2y >= 0, y <= x)
solve!(p)
```

* Indexing and Transpose
```
rows = 6
cols = 8
n = 2
X = Variable(rows, cols)
A = randn(rows, cols)
c = rand(1, n)
p = minimize(c * X[1:n, 5:5+n-1]' * c', X >= A)
solve!(p)
```

* Sum and Concatenation
```
x = Variable(4, 4)
y = Variable(4, 6)
p = maximize(sum(x) + sum(y), hcat(x, y) <= 2)
solve!(p)
```

* Minimum
```
x = Variable(10, 10)
y = Variable(10, 10)
a = rand(10, 10)
b = rand(10, 10)
p = maximize(minimum(min(x, y)), x <= a, y <= b)
solve!(p)
```

* Norm-Infinity
```
x = Variable(3)
p = minimize(norm(x, Inf), -2 <= x, x <= 1)
solve!(p)
```

## Solvers

By default, Convex.jl uses [ECOS](https://github.com/JuliaOpt/ECOS.jl) to solve SOCPs,
and [SCS](https://github.com/karanveerm/SCS.jl) to solve SDPs and exponential cone programs.
SCS currently works only on OSX, so SDPs and exponential cone programs
are only supported on OSX for now.
Any other solver in [JuliaOpt](https://github.com/JuliaOpt) may also be used, so long as it supports
the conic constraints used to express the problem. Currently,
other solvers will be most useful in solving linear programs (LPs)
and mixed integer linear programs (MILPs). Mosek AND Gurobi can be used to solve QPs.

For example, we can use GLPK to solve a MILP:
```
using GLPKMathProgInterface
solve!(p, GLPKSolverMIP())
```

You can set or see the current default solver by:
```
get_default_solver()
set_default_solver(GurobiSolver()) # or set_default_solver(ECOSSolver(verbose=0))
# Now Gurobi will be used by default as a solver
```

## Credits
Currently, Convex.jl is developed and maintained by:
- [Jenny Hong](http://www.stanford.edu/~jyunhong/)
- [Karanveer Mohan](http://www.stanford.edu/~kvmohan/)
- [Madeleine Udell](http://www.stanford.edu/~udell/)
- [David Zeng](http://www.stanford.edu/~dzeng0/)

The Convex.jl developers also thank:
- the [JuliaOpt](http://www.juliaopt.org/) team: [Iain Dunning](http://iaindunning.com/), [Joey Huchette](http://www.mit.edu/~huchette/) and [Miles Lubin](http://www.mit.edu/~mlubin/)
- [Stephen Boyd](http://www.stanford.edu/~boyd/), co-author of the book [Convex Optimization](http://www.stanford.edu/~boyd/books.html)
- [Steven Diamond](http://www.stanford.edu/~stevend2/), developer of [CVXPY](https://github.com/cvxgrp/cvxpy) and of a [DCP tutorial website](http://dcp.stanford.edu/) to teach disciplined convex programming.
- [Michael Grant](http://www.cvxr.com/bio), developer of [CVX](http://www.cvxr.com).

## Citing this package

If you use Convex.jl for published work,
we encourage you to cite the software using the following BibTeX citation:

    @article{convexjl,
     title = {Convex Optimization in {J}ulia},
     author ={Udell, Madeleine and Mohan, Karanveer and Zeng, David and Hong, Jenny and Diamond, Steven and Boyd, Stephen},
     year = {2014},
     journal = {SC14 Workshop on High Performance Technical Computing in Dynamic Languages},
     archivePrefix = "arXiv",
     eprint = {1410.4821},
     primaryClass = "math-oc",
    }
