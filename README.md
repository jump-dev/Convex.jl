# CVX.jl

<!--
[![Build Status](https://travis-ci.org/madeleineudell/CVX.jl.png)](https://travis-ci.org/madeleineudell/CVX.jl)
-->

CVX.jl is a julia package for disciplied convex programming.
This package is under active development; interfaces are not guaranteed to be stable, and bugs should be expected.
(Bug reports, of course, are welcome.)

CVX.JL allows you to express problems in simple, mathematical ways. All you need to worry about is the math, and CVX.jl will transfordm your problem into a standard form that is fed into a solver of your choice such as ECOS (and soon, SCS and solvers supported by MathProgBase.jl such as Gurobi, GLPK etc).

In its current state, CVX.jl can solve linear programs and very basic convex problems (max, min, norm_inf etc). For example, the following code solves a simple linear program:
```
x = Variable()
y = Variable()
constraints = [3x + y == 3, 4x + 3y >= 6, x + 2y <=3, x >=0, y >=0]
p = minimize(4x + y, constraints)
solve!(p)
p.optval # 3.6
x.value # 0.6
y.value # 1.2
```
# Prerequisites

CVX.jl requires
* [ECOS](http://github.com/ifa-ethz/ecos) >= 1.0.3

# Installation

To install, just git clone the repo. Work is still in progress and a clean, easy-to-use/install module will be available soon.

# Usage

A number of examples can be found in test/test.jl. Here are a few simple examples to start with:

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

