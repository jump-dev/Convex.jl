# Extended formulations and the DCP ruleset

Convex.jl works by transforming the problem (which possibly has nonsmooth,
nonlinear constructions like the nuclear norm, the log determinant, and so
forth—into) a linear optimization problem subject to conic constraints.

The transformed problem often involves adding auxiliary variables, and it is
called an "extended formulation," since the original problem has been extended
with additional variables.

Creating an extended formulation relies on the problem being modeled by
combining Convex.jl's "atoms" or primitives according to certain rules which
ensure convexity, called the
[disciplined convex programming (DCP) ruleset](http://cvxr.com/cvx/doc/dcp.html).
If these atoms are combined in a way that does not ensure convexity, the
extended formulations are often invalid.

## A valid formulation

As a simple example, consider the problem:
```@repl
using Convex, SCS
x = Variable();
model_min = minimize(abs(x), [x >= 1, x <= 2]);
solve!(model_min, SCS.Optimizer; silent_solver = true)
```

The optimum occurs at `x = 1`, but let us imagine we want to solve this problem
via Convex.jl using a linear programming (LP) solver.

Since `abs` is a nonlinear function, we need to reformulate the problem to pass
it to the LP solver. We do this by introducing an auxiliary variable `t` and
instead solving:
```@repl
using Convex, SCS
x = Variable();
t = Variable();
model_min_extended = minimize(t, [x >= 1, x <= 2, t >= x, t >= -x]);
solve!(model_min_extended, SCS.Optimizer; silent_solver = true)
```
That is, we add the constraints `t >= x` and `t >= -x`, and replace `abs(x)` by
`t`. Since we are minimizing over `t` and the smallest possible `t` satisfying
these constraints is the absolute value of `x`, we get the right answer. This
reformulation worked because we were minimizing `abs(x)`, and that is a valid
way to use the primitive `abs`.

## An invalid formulation

The reformulation of `abx(x)` works only if we are minimizing `t`.

Why? Well, let us consider the same reformulation for a maximization problem.
The original problem is:
```@repl
using Convex
x = Variable();
model_max = maximize(abs(x), [x >= 1, x <= 2])
```
This time, `problem is DCP` reports `false`. If we attempt to solve the problem,
an error is thrown:
```julia
julia> solve!(model_max, SCS.Optimizer; silent_solver = true)
┌ Warning: Problem not DCP compliant: objective is not DCP
└ @ Convex ~/.julia/dev/Convex/src/problems.jl:73
ERROR: DCPViolationError: Expression not DCP compliant. This either means that your problem is not convex, or that we could not prove it was convex using the rules of disciplined convex programming. For a list of supported operations, see https://jump.dev/Convex.jl/stable/operations/. For help writing your problem as a disciplined convex program, please post a reproducible example on https://jump.dev/forum.
Stacktrace:
[...]
```

The error is thrown because, if we do the same reformulation as before, we
arrive at the problem:
```@repl
using Convex
x = Variable();
t = Variable();
model_max_extended = maximize(t, [x >= 1, x <= 2, t >= x, t >= -x]);
solve!(model_max_extended, SCS.Optimizer; silent_solver = true)
```
whose solution is unbounded.

In other words, we can get the wrong answer by using the extended reformulation,
because the extended formulation was only valid for a minimization problem.

Convex.jl always creates the extended reformulation, but because they are only
guaranteed to be valid when the DCP ruleset is followed, Convex.jl will
programmatically check the whether or not these DCP rules were satisfied and
error if they were not.
