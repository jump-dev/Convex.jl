# Convex.jl - Convex Optimization in Julia

Convex.jl is a Julia package for [Disciplined Convex Programming](http://dcp.stanford.edu/) (DCP).

Convex.jl makes it easy to describe optimization problems in a natural,
mathematical syntax, and to solve those problems using a variety of different
(commercial and open-source) solvers.

Convex.jl can be used to solve:

 * linear programs
 * mixed-integer linear programs and mixed-integer second-order cone programs
 * DCP-compliant convex programs including
   * second-order cone programs (SOCP)
   * exponential cone programs
   * semidefinite programs (SDP)

## Resources for getting started

There are a few ways to get started with Convex:

 * Read the [Installation](@ref) guide
 * Read the introductory tutorial [Quick Tutorial](@ref)
 * Read the list of [Supported Operations](@ref)
 * Browse some of our examples

!!! tip
    Need help? Join the [community forum](https://jump.dev/forum)
    to search for answers to commonly asked questions.

    Before asking a question, make sure to read the post [make it easier to help you](https://discourse.julialang.org/t/psa-make-it-easier-to-help-you/14757),
    which contains a number of tips on how to ask a good question.

## How the documentation is structured

Having a high-level overview of how this documentation is structured will help
you know where to look for certain things.

* **Examples** contain worked examples of solving problems with Convex. Start
  here if you are new to Convex, or you have a particular problem class you want
  to model.

* The **Manual** contains short code-snippets that explain how to achieve
  specific tasks in Convex. Look here if you want to know how to achieve a
  particular task.

* The **Developer docs** section contains information for people contributing to
  Convex development. Don't worry about this section if you are using Convex to
  formulate and solve problems as a user.

## Extended formulations and the DCP ruleset

Convex.jl works by transforming the problem (which possibly has nonsmooth,
nonlinear constructions like the nuclear norm, the log determinant, and so
forth—into) a linear optimization problem subject to conic constraints. This
reformulation often involves adding auxiliary variables, and is called an
"extended formulation," since the original problem has been extended with
additional variables. These formulations rely on the problem being modeled by
combining Convex.jl's "atoms" or primitives according to certain rules which
ensure convexity, called the
[disciplined convex programming (DCP) ruleset](http://cvxr.com/cvx/doc/dcp.html).
If these atoms are combined in a way that does not ensure convexity, the
extended formulations are often invalid. As a simple example, consider the problem

```julia
model = minimize(abs(x), x >= 1, x <= 2)
```

The optimum occurs at `x=1`, but let us imagine we want to solve this problem
via Convex.jl using a linear programming (LP) solver.

Since `abs` is a nonlinear function, we need to reformulate the problem to pass
it to the LP solver. We do this by introducing an auxiliary variable `t` and
instead solving:
```julia
model = minimize(t, x >= 1, x <= 2, t >= x, t >= -x)
```
That is, we add the constraints `t >= x` and `t >= -x`, and replace `abs(x)` by
`t`. Since we are minimizing over `t` and the smallest possible `t` satisfying
these constraints is the absolute value of `x`, we get the right answer. This
reformulation worked because we were minimizing `abs(x)`, and that is a valid
way to use the primitive `abs`.

If we were maximizing `abs`, Convex.jl would error with

> Problem not DCP compliant: objective is not DCP

Why? Well, let us consider the same reformulation for a maximization problem.
The original problem is:
```julia
model = maximize(abs(x), x >= 1, x <= 2)
```
and the maximum of 2, obtained at `x = 2`. If we do the same reformulation as
above, however, we arrive at the problem:
```julia
maximize(t, x >= 1, x <= 2, t >= x, t >= -x)
```
whose solution is infinity.

In other words, we got the wrong answer by using the reformulation, since the
extended formulation was only valid for a minimization problem. Convex.jl always
performs these reformulations, but they are only guaranteed to be valid when the
DCP ruleset is followed. Therefore, Convex.jl programmatically checks the
whether or not these rules were satisfied and errors if they were not.
