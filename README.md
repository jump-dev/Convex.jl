# CVX.jl

<!--
[![Build Status](https://travis-ci.org/madeleineudell/CVX.jl.png)](https://travis-ci.org/madeleineudell/CVX.jl)
-->

CVX.jl is a julia package for [Disciplined Convex Programming](http://dcp.stanford.edu/). This package is under active development; interfaces are not guaranteed to be stable, and bugs should be expected. We'd love bug reports and feature requests!

## Table of Contents
- [Introduction](#user-content-introduction)
- [Supported Operations](#user-content-supported-operations)
- [Installation](#user-content-installation)
- [Basic Types](#user-content-basic-types)
  - [Variables](#user-content-variables)
  - [Constants](#user-content-constants)
  - [Expressions](#user-content-expressions)
  - [Constraints](#user-content-constraints)
  - [Objective](#user-content-objective)
  - [Problem](#user-content-problem)
- [Examples](#user-content-examples)
- [Credits](#user-content-credits)


## Introduction
CVX.jl allows you to express problems in simple, mathematical ways. All you need to worry about is the math, and CVX.jl will transform your problem into a standard form that is fed into a solver of your choice such as ECOS (and soon, SCS). Here's a quick example of code that solves a least-squares problem with inequality constraints:

```
# Generate problem data
m = 4
n = 5
A = randn(m, n)
b = randn(m, 1)

# Create a variable of size n x 1. 
# Matrix variables are also supported
x = Variable(n)

# The problem is to minimize ||Ax + b||^2 subject to x >= 1
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

# The dual values are stored with the corresponding constraints
problem.constraints[1].dual_value
```

## Supported Operations
In its current state, CVX.jl supports affine constraints and second-order cone (SOC) constraints. In most cases these have been suitably overloaded to work seamlessly for scalars, vectors and matrices. A list of operations that can be performed are listed below:

- Affine
 - addition, subtraction, multiplication, division: `+, -, /, *`
 - indexing into vectors and matrices: `x[1:4, 2:3]`
 - k-th diagonal of a matrix: `diag(x, k)`
 - transpose: `x'`
 - dot product: `x' * y` or `dot(x, y)`
 - reshape, vec: `reshape(x, 2, 3)` or `vec(x)`
 - min, max element of a vector or matrix: `max(x)`
 - horizontal and vertical stacking: `hcat(x, y); vertcat(x, y)`
- Elementwise
 - elementwise min, max between vectors or matrices: `max(x, y)`
 - pos, neg where pos is implemented as max(x, 0) and neg is implemented as -max(-x, 0): `pos(x)`
 - inverse pos (1./pos(x)): `inv_pos(x)`
 - square root: `sqrt(x)`
 - square: `square(x)`
 - square_pos (square(pos(x)): `square_pos(x)`
 - absolute value: `abs(x)`
- SOC/ Other supported constraints
 - geometric mean: `geo_mean(x, y)`
 - norm (norm_1, norm_2, norm_inf): `norm(x, 1); norm(x, 2); norm(x, Inf)`
 - quadratic form: `quad_form(P, x)`
 - quadratic over linear: `quad_over_lin(x, y)`
 - l2-norm squared: `sum_squares(x)`

In addition to these operations, when there are operations between vectors/matrices and scalars, the scalars are promoted. For example, we can do `max(x, 0)` where x is a vector variable but 0 is a scalar.


## Installation
```
Pkg.clone("https://github.com/cvxgrp/CVX.jl.git")
Pkg.clone("git@github.com:karanveerm/ECOS.jl.git")
Pkg.build("ECOS")  
```
You might have to restart Julia.
This does not work on Windows yet but should work on Mac/ Linux. Please file an issue in case you run into problems during installation. We'll be glad to help!

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
CVX.jl allows the values of the expressions to be evaluated directly. 
```
x = Variable()
y = Variable()
z = Variable()
expr = x + y + z
problem = minimize(expr, [ x >= 1, y >= x, 4 * z >= y]
solve!(problem)

# Once the problem is solved, we can call evaluate() on expr:
expr.evaluate()
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

After the problem is solved, the status can be checked by `problem.status`, which can be `solved`, `primal infeasible, `dual infeasible`, `max iterations reached` or `numerical problems in solver`. If the status is `solved`, `problem.optval` will have the optimum value of the problem. Each variable has a `value` that can be used to access the variables optimum value. The optimum value of expressions can also be found by calling the `evaluate()` function of the expression as follows: `expr.evaluate()`. The dual values are stored with the respective constraints and can be accessed as `problem.constraints[idx].dual_value`.

## Examples
A number of very simple can be found in test/test.jl. More sophisticated examples, along with plots can be found in the examples/ directory.
Here are a few simple examples to start with:

* Dot Product
```
x = Variable(2)
A = 1.5 * eye(2)
p = minimize(dot([2.0; 2.0], x), [A * x >= [1.1; 1.1]])
solve!(p)
println(p.optval)
println(x.value)
```

* Matrix Variables
```
X = Variable(2, 2)
c = ones(2, 1)
p = minimize(c' * X * c, [X >= ones(2, 2)])
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
p = minimize(objective, [x >= 3, 2y >= 0, y <= x])
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
p = maximize(sum(x) + sum(y), [hcat(x, y) <= 2])
solve!(p)
```

* Minimum
```
x = Variable(10, 10)
y = Variable(10, 10)
a = rand(10, 10)
b = rand(10, 10)
p = maximize(min(min(x, y)), [x <= a, y <= b])
solve!(p)
```

* Norm-Infinity
```
x = Variable(3)
p = minimize(norm(x, Inf), [-2 <= x, x <= 1])
solve!(p)
```

## Credits
Currently, CVX.jl is developed and maintained by:
- [Jenny Hong](http://www.stanford.edu/~jyunhong/)
- [Karanveer Mohan](http://www.stanford.edu/~kvmohan/)
- [Madeleine Udell](http://www.stanford.edu/~udell/)
- [David Zeng](http://www.stanford.edu/~dzeng0/)

In addition to development, we'd like to give a huge thanks to:
- [Stephen Boyd](http://www.stanford.edu/~boyd/): Professor of Electrical Engineering, Stanford University. He is also the co-author of the book [Convex Optimization](http://www.stanford.edu/~boyd/books.html). We thank Professor Boyd for his continuous input and support.
- [Steven Diamond](http://www.stanford.edu/~stevend2/): CVX.jl started out as a wrapper around Steven Diamond's [CVXPY](https://github.com/cvxgrp/cvxpy) and its design has been inspired from CVXPY. We greatly appreciate Steven Diamond's experienced and continual guidance. In addition, Steven Diamond also wrote http://dcp.stanford.edu/ to teach disciplined convex programming, a useful resource for CVX.jl users.
