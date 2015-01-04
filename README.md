# Convex.jl

[![Build Status](https://travis-ci.org/cvxgrp/Convex.jl.svg?branch=master)](https://travis-ci.org/cvxgrp/Convex.jl) [![Coverage Status](https://img.shields.io/coveralls/cvxgrp/Convex.jl.svg)](https://coveralls.io/r/cvxgrp/Convex.jl)

Convex.jl is a julia package for [Disciplined Convex Programming](http://dcp.stanford.edu/). Convex.jl can solve linear programs, mixed-integer linear programs, and dcp-compliant convex programs using a variety of solvers, including [Mosek](https://github.com/JuliaOpt/Mosek.jl), [Gurobi](https://github.com/JuliaOpt/gurobi.jl), [ECOS](https://github.com/JuliaOpt/ECOS.jl), [SCS](https://github.com/karanveerm/SCS.jl), [GLPK](https://github.com/JuliaOpt/GLPK.jl), through the [MathProgBase](http://mathprogbasejl.readthedocs.org/en/latest/) interface.

Detailed documentation, installation instructions and examples for Convex.jl can be found [here](http://convexjl.readthedocs.org/). Below, we give some basic usage guide.


## Installation
Convex.jl can be installed on Julia version 0.3 or above
```
Pkg.add("Convex")
```

## Quick Example
Here's a quick example of code that solves a least-squares problem with inequality constraints
```
# Let us first make the Convex.jl module available
using Convex

# Generate random problem data
m = 4;  n = 5
A = randn(m, n); b = randn(m, 1)

# Create a (column vector) variable of size n x 1.
x = Variable(n) # or x = Variable(n, 1)

# The problem is to minimize ||Ax - b||^2 subject to x >= 0
# This can be done by: minimize(objective, constraints)
problem = minimize(sum_squares(A * x + b), [x >= 0])

# You can add more constraints at any time by
problem.constraints += [x <= 1, 0.5 <= 2 * x]

# Solve the problem by calling solve!
solve!(problem)

# Check the status of the problem
problem.status # :Optimal, :Infeasible, :Unbounded etc.

# Get the optimum value
problem.optval
```

Note that Convex.jl was previously called CVX.jl. If you're running into bugs or have feature requests, please use the [Github Issue Tracker](https://github.com/cvxgrp/Convex.jl/issues>). For usage questions, please contact us via the [JuliaOpt mailing list](https://groups.google.com/forum/#!forum/julia-opt>).


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

