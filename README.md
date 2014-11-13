# Convex.jl

<!--
[![Build Status](https://travis-ci.org/cvxgrp/Convex.jl.png)](https://travis-ci.org/cvxgrp/Convex.jl)
-->

Convex.jl is a julia package for [Disciplined Convex Programming](http://dcp.stanford.edu/). Note that Convex.jl was previously called CVX.jl. This package is under active development; interfaces are not guaranteed to be stable, and bugs should be expected. Nevertheless, we try to fix problems that come up as swiftly as we can. We'd love bug reports and feature requests!

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
- [Credits](#user-content-credits)


## Introduction
Convex.jl allows you to express problems in simple, mathematical ways. All you need to worry about is the math, and Convex.jl will transform your problem into a standard form that is fed into a solver of your choice such as ECOS (and soon, SCS). Here's a quick example of code that solves a least-squares problem with inequality constraints:

```
using Convex

# Generate problem data
m = 4
n = 5
A = randn(m, n)
b = randn(m, 1)

# Create a variable of size n x 1.
# Matrix variables are also supported
x = Variable(n)

# The problem is to minimize ||Ax + b||^2 subject to x >= 0
problem = minimize(sum_squares(A * x + b), [0 <= x])

# Alternatively, we can add constraints at any time in the following way:
problem.constraints += [x <= 1, 0.5 <= 2*x]

# Solve the problem by calling solve!
solve!(problem)

# Status (solved, infeasible etc.)
problem.status

# Optimum value
problem.optval

# Optimal value of x
x.value
```

## Installation
```
Pkg.clone("https://github.com/cvxgrp/Convex.jl.git")
Pkg.clone("https://github.com/karanveerm/ECOS.jl.git")
Pkg.build("ECOS")
```
If you're on OSX, then SCS.jl should also work. This can be used to solve problems involving exponential and semi-definite constraints.
NOTE: SCS.jl ONLY works on OSX.
```
Pkg.clone("https://github.com/karanveerm/SCS.jl.git")
Pkg.build("SCS")
```
You might have to restart Julia.
Please file an issue in case you run into problems during installation. We'll be glad to help!

## Basic Types

### Variables
Variables represent the quantities that we want to find by solving the problem.
```
# Scalar variable
x = Variable()
y = Variable(1)

# Column-vector variable
x = Variable(5)

# Matrix variable
x = Variable(4, 6)
```

### Constants
Use constants the way you would usually use them in Julia. They should work just fine. In case there are any problems, please report an issue and we will fix it as soon as we can.

### Expressions
Performing operations on Variables results in Expressions. As the name suggests, these are mathematical expressions. These expressions can be combined with other expressions and so on. Expressions that are created must be DCP-compliant (follow the rules of Disciplined Convex Programming (DCP)). More information on DCP can be found here: http://dcp.stanford.edu/.
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
Constants, Variables and Expressions,  combined with <, >, <=, >= and == create Constraints.
```
x = Variable(5, 5)
# Equality constraint
x == 0
# Greater than/ equal to constraint
x >= 1
```
If a constraint is specified as strictly greater than or less than, it will be converted to a >= or <= constraint respectively. If a strict inequality is needed, one way to do it would be as follows:
```
x = Variable()
slack = Variable()
constraints = [x + slack >= 0, slack >= 1e-4]
```

### Objective
The objective of the problem is a scalar expression to be maximized or minimized by using `maximize` or `minimize` respectively. Feasibility problems are also allowed by either giving a constant as the expression, or using `problem = satisfy(constraints)`.

### Problem
A problem is an objective and a list of constraints. These are constructed using the form
```problem = minimize(objective, constraints)```
or
```problem = maximize(objective, constraints)```
The constraints can be added at any time before the problem is solved.
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

After the problem is solved, the status can be checked by `problem.status`, which can be `:Optimal`, `:Infeasible`, `:Unbounded`, `:Indeterminate` or `:Error`. If the status is `:Optimal`, `problem.optval` will have the optimum value of the problem. Each variable has a `value` that can be used to access the variables optimum value. The optimum value of expressions can also be found by calling the `evaluate()` function of the expression as follows: `evaluate(expr)`. <!--The dual values are stored with the respective constraints and can be accessed as `problem.constraints[idx].dual_value`.-->


## Supported operations
Convex.jl currently supports the following operations. Except where explicitly noted below, they work seamlessly for scalars, vectors and matrices. These atomic operations ("atoms") may be composed according to the [DCP](dcp.stanford.edu) composition rules to form new convex, concave, or affine expressions.

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

### Vector and Matrix atoms

These atoms take vector or matrix arguments and return scalar expressions.

atom | description | vexity | slope | implicit constraint
-----|--------|--------|------- |-------
`norm(x, p)` | $(\sum x_i^p)^{1/p}$ | convex | increasing on $x \ge 0$, decreasing on $x \le 0$ | `p = 1, 2, Inf`
`vecnorm(x, p)` | $(\sum x_{ij}^p)^{1/p}$ | convex | increasing on $x \ge 0$, decreasing on $x \le 0$ | `p = 1, 2, Inf`
`quad_form(x,P)` | $x^T P x$ | convex in $X$, affine in $P$ | increasing on $x \ge 0$, decreasing on $x \le 0$, increasing in $P$ | either $x$ or $P$ must be constant
`quad_over_lin(x, y)` | $x^T x/y$ | convex | increasing on $x \ge 0$, decreasing on $x \le 0$, decreasing in $y$ | $y > 0$
`sum_squares(x)` | $\sum x_i^2$ | convex | increasing on $x \ge 0$, decreasing on $x \le 0$ | none

### Promotion
When an atom or constraint is applied to a scalar and a higher dimensional variable, the scalars are promoted. For example, we can do `max(x, 0)` gives an expression with the shape of `x` whose elements are the maximum of the corresponding element of `x` and `0`.

## Examples
A number of very simple can be found in test/test.jl. More sophisticated examples, along with plots can be found in the examples/ directory.
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

## Credits
Currently, Convex.jl is developed and maintained by:
- [Jenny Hong](http://www.stanford.edu/~jyunhong/)
- [Karanveer Mohan](http://www.stanford.edu/~kvmohan/)
- [Madeleine Udell](http://www.stanford.edu/~udell/)
- [David Zeng](http://www.stanford.edu/~dzeng0/)

The Convex.jl developers also thank:
- [Stephen Boyd](http://www.stanford.edu/~boyd/): Professor of Electrical Engineering, Stanford University. He is also the co-author of the book [Convex Optimization](http://www.stanford.edu/~boyd/books.html). We thank Professor Boyd for his continuous input and support.
- [Steven Diamond](http://www.stanford.edu/~stevend2/): many aspects of the design of Convex.jl were inspired Steven Diamond's [CVXPY](https://github.com/cvxgrp/cvxpy). We greatly appreciate Steven Diamond's experienced and continual guidance. In addition, Steven Diamond also wrote a [DCP tutorial website](http://dcp.stanford.edu/) to teach disciplined convex programming, a useful resource for Convex.jl users.

## Citing this package

If you use Convex.jl for published work, 
we encourage you to cite the software.

Use the following BibTeX citation:

    @article{udell2014,
        title = {Convex Optimization in Julia},
        author ={Udell, Madeleine and Mohan, Karanveer and Zeng, David and Hong, Jenny and Diamond, Steven and Boyd, Stephen},
        year = {2014},
        archivePrefix = "arXiv",
        eprint = {1410.4821},
        primaryClass = "math-oc",
        journal={arXiv preprint arXiv:1410.4821},
    }
