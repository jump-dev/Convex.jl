# Convex.jl

[![Build Status](https://travis-ci.org/JuliaOpt/Convex.jl.svg?branch=master)](https://travis-ci.org/JuliaOpt/Convex.jl) [![Coverage Status](https://img.shields.io/coveralls/JuliaOpt/Convex.jl.svg)](https://coveralls.io/r/JuliaOpt/Convex.jl)

Convex.jl is a julia package for [Disciplined Convex Programming](http://dcp.stanford.edu/). Convex.jl can solve linear programs, mixed-integer linear programs, and dcp-compliant convex programs using a variety of solvers, including [Mosek](https://github.com/JuliaOpt/Mosek.jl), [Gurobi](https://github.com/JuliaOpt/gurobi.jl), [ECOS](https://github.com/JuliaOpt/ECOS.jl), [SCS](https://github.com/karanveerm/SCS.jl), [GLPK](https://github.com/JuliaOpt/GLPK.jl), through the [MathProgBase](http://mathprogbasejl.readthedocs.org/en/latest/) interface.

*Detailed documentation, installation instructions and examples for Convex.jl can be found [here](http://convexjl.readthedocs.org/).* 
Below, we give a brief introduction to basic usage.


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
x = Variable(n)

# The problem is to minimize ||Ax - b||^2 subject to x >= 0
# This can be done by: minimize(objective, constraints)
problem = minimize(sum_squares(A * x + b), [x >= 0])

# Solve the problem by calling solve!
solve!(problem)

# Check the status of the problem
problem.status # :Optimal, :Infeasible, :Unbounded etc.

# Get the optimal value
problem.optval
```

Note that Convex.jl was previously called CVX.jl. If you're running into bugs or have feature requests, please use the [Github Issue Tracker](https://github.com/JuliaOpt/Convex.jl/issues>). For usage questions, please contact us via the [JuliaOpt mailing list](https://groups.google.com/forum/#!forum/julia-opt>).


## Examples
A number of examples can be found [here](http://convexjl.readthedocs.org/en/latest/examples.html). 
The [basic usage notebook](http://nbviewer.ipython.org/github/JuliaOpt/Convex.jl/blob/master/examples/basic_usage.ipynb) gives a simple tutorial on problems that can be solved using Convex.jl


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
