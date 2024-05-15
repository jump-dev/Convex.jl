# FAQ

## Where can I get help?

For usage questions, please start a new post on the
[Julia Discourse](https://discourse.julialang.org/c/domain/opt).

If you have a reproducible example of a bug or if you have a feature request,
please open a [GitHub issue](https://github.com/jump-dev/Convex.jl/issues/new).

## How does Convex.jl differ from JuMP?

Convex.jl and JuMP are both modelling languages for mathematical programming
embedded in Julia, and both interface with solvers via
[MathOptInterface](https://github.com/jump-dev/MathOptInterface.jl).

Convex.jl converts problems to a standard conic form. This approach requires
(and certifies) that the problem is convex and DCP compliant, and guarantees
global optimality of the resulting solution (if the solver succeeds. For some
models, the solver may experience numerical difficulty).

JuMP allows nonlinear programming through an interface that learns about
functions via their derivatives. This approach is more flexible (for example,
you can optimize non-convex functions), but can't guarantee global optimality if
your function is not convex, or warn you if you've entered a non-convex
formulation.

For linear programming, the difference is more stylistic. JuMP's syntax is
scalar-based and similar to AMPL and GAMS making it easy and fast to create
constraints by indexing and summation (like `sum(x[i] for i in 1:n)`).

Convex.jl allows (and prioritizes) linear algebraic and functional constructions
(like `max(x, y) <= A * z`); indexing and summation are also supported in Convex.jl,
but are somewhat slower than in JuMP.

JuMP also lets you efficiently solve a sequence of problems when new constraints
are added or when coefficients are modified, whereas Convex.jl parses the
problem again whenever the [solve!](@ref) method is called.

## Where can I learn more about Convex Optimization?

See the freely available book [Convex Optimization](http://web.stanford.edu/~boyd/cvxbook/)
by Boyd and Vandenberghe for general background on convex optimization.

For help understanding the rules of Disciplined Convex Programming, see the
[DCP tutorial website](http://dcp.stanford.edu/).

## Are there similar packages available in other languages?

See [CVXPY](http://www.cvxpy.org) in Python and [CVX](http://cvxr.com/) in
Matlab.

## How does Convex.jl work?

For a detailed discussion of how Convex.jl works, see [our paper](http://www.arxiv.org/abs/1410.4821).

## How do I cite this package?

If you use Convex.jl for published work, we encourage you to cite the software
using the following BibTeX citation:

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
