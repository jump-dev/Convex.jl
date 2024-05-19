# Performance tips

There are 3 phases of building and solving a Convex.jl problem:

1. Building the expression tree. This is everything that happens before `solve!` is called, e.g. creating variables, problems, constraints, and other expressions.
2. Formulating the problem in MathOptInterface. This happens automatically when `solve!` is called, and (as long as `solve!` is called with `silent=false`) emits a log message when completed, containing the amount of time it took and the total amount of memory allocations.
3. Solving the problem. This happens in `solve!` after the problem is formulated. When `silent=false`, typically solvers will print information about the problem and their progress when solving it.

If you are experiencing performance issues, it is important to isolate which of these three steps is taking the most time. Typically, step (1) should be the fastest, then step (2), then step (3).

* If building the expression tree (step 1 above) is the slowest part, there is likely a performance bug in Convex.jl. Please [file an issue](https://github.com/jump-dev/Convex.jl/issues/new) with a reproducible minimal example; it may be able to be quickly fixed.
* If formulating the problem is much slower than solving it, there may be a performance issue in Convex.jl. See [Faster problem formulation: Avoid scalar indexing](@ref) below for one tip to speed things up.
* If solving the problem is the slowest part, you may still be able to improve performance. See [Speeding up solve times: dualization](@ref) and [Speeding up solve times: conic-friendly formulations](@ref) below. Also consider trying other solvers. For large problems, first-order solvers like SCS and COSMO may be more performant (albeit with higher convergence tolerances) than second-order solvers like Clarabel and Hypatia.

If you're struggling to improve performance, feel free to ask for help in the [community forum](https://jump.dev/forum). Before asking a question, make sure to read the post [make it easier to help you](https://discourse.julialang.org/t/psa-make-it-easier-to-help-you/14757), which contains a number of tips on how to ask a good question.

## Faster problem formulation: Avoid scalar indexing

When manipulating Convex.jl expressions, try to use a "vectorized" style, and avoid loops and scalar indexing. For example,

```julia
function my_dot(x, y) # avoid this kind of code for Convex.jl expressions!
    s = 0.0
    for i in eachindex(x, y)
        s += x[i]*y[i]
    end
    return s
end
```

is usually perfectly reasonable Julia code (although perhaps `s` should be initialized as `s=zero(promote_type(eltype(x), eltype(y))))` to be more generic), since loops are typically fast in Julia.

However, when using Convex.jl expressions like `x = Variable(100)` this will create 100 separate `IndexAtom`s, each of which has some overhead! That is because Convex maintains a lazy expression tree representing the problem. Instead, `sum(x .* y)` is much more efficient for Convex.jl problem formulation. In this case of course, one can just use `dot` directly.

Indexing style typically should have no effect on solve times (step (3) above), only on constructing the expression tree and formulating the problem.

## Speeding up solve times: dualization

Depending on the solver and the problem, it can sometimes speed things up quite a lot to pass the solver the dual problem to solve rather than the primal problem. This can be accomplished easily with [Dualization.jl](https://github.com/jump-dev/Dualization.jl), simply by passing `dual_optimizer(optimizer)` instead of `optimizer` to `solve!`. See [Dualization](@ref) for an example.

## Speeding up solve times: conic-friendly formulations

Convex.jl is built around conic programming, and works with conic solvers. This is a different optimization methodology than you may be used to, and works by formulating problems in terms of affine objective functions and affine-function-in-cone constraints for a variety of convex cones. Convex.jl reformulates all problems to this form. Specifically, the e.g. objective function you give Convex will not actually be executed line-by-line; instead, it will be reformulated to a possibly very-different looking form for the conic solver to handle. This has important consequences for performance.

For example, it is typically easier for conic solvers to minimize $\|A*x\|$ than $\|A*x\|^2$, since the former is reformulated as a second order cone constraint, and the latter as the same second order cone constraint, plus a rotated second order cone constraint (to model the final squaring). This may be unintuitive, since for non-conic solvers, often $\|A*x\|^2$ is easier to solve (as it avoids a square-root in the typical formulation).

Convex.jl provides an abstraction by automatically performing these reformulations, but as often is the case, performance leaks through. To have the best performance, you may will likely need to learn more details about the various cones supported by your solver, by [MathOptInterface](https://jump.dev/MathOptInterface.jl/stable/reference/standard_form/#Vector-sets), and how Convex.jl reformulates your problem.
