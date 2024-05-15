## Extended formulations and the DCP ruleset

Convex.jl works by transforming the problem (which possibly has nonsmooth,
nonlinear constructions like the nuclear norm, the log determinant, and so
forth—into) a linear optimization problem subject to conic constraints. This
reformulation often involves adding auxiliary variables, and is called an
"extended formulation," since the original problem has been extended with
additional variables.

Creating an extended formulation relies on the problem being modeled by
combining Convex.jl's "atoms" or primitives according to certain rules which
ensure convexity, called the
[disciplined convex programming (DCP) ruleset](http://cvxr.com/cvx/doc/dcp.html).
If these atoms are combined in a way that does not ensure convexity, the
extended formulations are often invalid.

As a simple example, consider the problem:
```@repl
using Convex
x = Variable();
model = minimize(abs(x), [x >= 1, x <= 2])
```

The optimum occurs at `x = 1`, but let us imagine we want to solve this problem
via Convex.jl using a linear programming (LP) solver.

Since `abs` is a nonlinear function, we need to reformulate the problem to pass
it to the LP solver. We do this by introducing an auxiliary variable `t` and
instead solving:
```@repl
using Convex
x = Variable();
t = Variable();
model = minimize(t, [x >= 1, x <= 2, t >= x, t >= -x])
```
That is, we add the constraints `t >= x` and `t >= -x`, and replace `abs(x)` by
`t`. Since we are minimizing over `t` and the smallest possible `t` satisfying
these constraints is the absolute value of `x`, we get the right answer. This
reformulation worked because we were minimizing `abs(x)`, and that is a valid
way to use the primitive `abs`.

This reformulation works only if we are minimizing `t`.

Why? Well, let us consider the same reformulation for a maximization problem.
The original problem is:
```@repl
using Convex
x = Variable();
model = maximize(abs(x), [x >= 1, x <= 2])
```
and the maximum of 2, obtained at `x = 2`. If we do the same reformulation as
above, however, we arrive at the problem:
```@repl
using Convex
x = Variable();
t = Variable();
maximize(t, [x >= 1, x <= 2, t >= x, t >= -x])
```
whose solution is infinity.

In other words, we got the wrong answer by using the extended reformulation,
because the extended formulation was only valid for a minimization problem.

Convex.jl always performs these reformulations, but they are only guaranteed to
be valid when the DCP ruleset is followed. Therefore, Convex.jl programmatically
checks the whether or not these rules were satisfied and errors if they were not.
