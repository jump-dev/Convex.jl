Convex.jl - Convex Optimization in Julia
========================================

Convex.jl is a Julia package for [Disciplined Convex
Programming](http://dcp.stanford.edu/) (DCP). Convex.jl makes it easy to
describe optimization problems in a natural, mathematical syntax, and to solve
those problems using a variety of different (commercial and open-source)
solvers. Convex.jl can solve

-   linear programs
-   mixed-integer linear programs and mixed-integer second-order cone programs
-   dcp-compliant convex programs including
    -   second-order cone programs (SOCP)
    -   exponential cone programs
    -   semidefinite programs (SDP)

Convex.jl supports many solvers, including
[COSMO](https://github.com/oxfordcontrol/COSMO.jl),
[Mosek](https://github.com/JuliaOpt/Mosek.jl),
[Gurobi](https://github.com/JuliaOpt/gurobi.jl),
[ECOS](https://github.com/JuliaOpt/ECOS.jl),
[SCS](https://github.com/karanveerm/SCS.jl) and
[GLPK](https://github.com/JuliaOpt/GLPK.jl), through
[MathOptInterface](https://github.com/JuliaOpt/MathOptInterface.jl).

Note that Convex.jl was previously called CVX.jl. This package is under active
development; we welcome bug reports and feature requests. For usage questions,
please contact us via the [Julia
Discourse](https://discourse.julialang.org/c/domain/opt).

## Extended formulations and the DCP ruleset

Convex.jl works by transforming the problem—which possibly has nonsmooth,
nonlinear constructions like the nuclear norm, the log determinant, and so
forth—into a linear optimization problem subject to conic constraints. This
reformulation often involves adding auxiliary variables, and is called an
"extended formulation", since the original problem has been extended with
additional variables. These formulations rely on the problem being modelled by
combining Convex.jl's "atoms" or primitives according to certain rules which
ensure convexity, called the [disciplined convex programming (DCP)
ruleset](http://cvxr.com/cvx/doc/dcp.html). If these atoms are combined in a way
that does not ensure convexity, the extended formulations are often invalid. As
a simple example, consider the problem

```julia
minimize( abs(x), x >= 1, x <= 2)
```

Obviously, the optimum occurs at `x=1`, but let us imagine we want to solve this
problem via Convex.jl using a linear programming (LP) solver. Since `abs` is a
nonlinear function, we need to reformulate the problem to pass it to the LP
solver. We do this by introducing an auxiliary variable `t` and instead solving

```julia
minimize(t, x >= 1, x <= 2, t >= x, t >= -x)
```

That is, we add the constraints `t >= x` and `t >= -x`, and replace `abs(x)` by
`t`. Since we are minimizing over `t` and the smallest possible `t` satisfying
these constraints is the absolute value of `x`, we get the right answer. That
is, this reformulation worked because we were minimizing `abs(x)`, and that is a
valid way to use the primitive `abs`.

If we were maximizing `abs`, Convex.jl would print

> Warning: Problem not DCP compliant: objective is not DCP

Why? Well, let us consider the same reformulation for a maximization problem.
The original problem is now

```julia
maximize( abs(x), x >= 1, x <= 2)
```

and trivially the optimum is 2, obtained at `x=2`. If we do the same
replacements as above, however, we arrive at the problem

```julia
maximize(t, x >= 1, x <= 2, t >= x, t >= -x)
```

whose solution is infinity. In other words, we got the wrong answer by using the
reformulation, since the extended formulation was only valid for a minimization
problem. Convex.jl always performs these reformulations, but they are only
guaranteed to be valid when the DCP ruleset is followed. Therefore, Convex.jl
programatically checks the whether or not these rules were satisfied and warns
if they were not. One should not take these DCP warnings lightly!
