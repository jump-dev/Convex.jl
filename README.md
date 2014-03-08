# CVX

[![Build Status](https://travis-ci.org/madeleineudell/CVX.jl.png)](https://travis-ci.org/madeleineudell/CVX.jl)

CVX.jl is a julia package for disciplied convex programming.
This package is under active development; interfaces are not guaranteed to be stable, and bugs should be expected.
(Bug reports, of course, are welcome.)

# Installation

To install, just open a Julia prompt and call

    Pkg.clone("git@github.com:madeleineudell/CVX.jl.git")

You'll also need to install the python module cvxpy and its dependencies, cvxopt and ecos.

    easy_install cvxopt; easy_install ecos; easy_install cvxpy

# Usage

The interface for this package is not yet stable, but we give a few usage examples to start with.

* Norm minimization
```
y = Variable(2,2);
p = minimize(norm(y) - 1,
    y >= 1)
println(p)
solve!(p)
println(p.optval)
println(y.value)
```

* Quadratic programming
```
y = Variable(2);
A = randn(4,2);
b = randn(4);
p = minimize(sum(square(A*y-b)),
    y >= 0)
println(p)
solve!(p)
println(p.optval)
println(y.value)
```

* Nonnegative least squares, using an augmented lagrangian to illustrate functional notation
```
n = 4; m=5; rho=1
x = Variable(n); z = Variable(n); y = zeros(n);
D = randn(m,n);
h = randn(m);
A = eye(n);
B = -eye(n);
c = zeros(n);
f(x) = sum(square(D*x - h))
L(x,y,z,rho) = f(x) + rho/2*sum(square(A*x+B*z-c+y))
solve!(minimize(L(x,y,z),z>=0,A*x+B*z==c))
```
