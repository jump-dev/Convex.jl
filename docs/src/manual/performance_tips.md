# Performance tips

There are 3 phases of building and solving a Convex.jl problem:

1. Building the expression tree. This is everything that happens before `solve!` is called, e.g. creating variables, problems, constraints, and other expressions.
2. Formulating the problem in MathOptInterface. This happens automatically when `solve!` is called, and (as long as `solve!` is called with `silent=false`) emits a log message when completed, containing the amount of time it took and the total amount of memory allocations.
3. Solving the problem. This happens in `solve!` after the problem is formulated. When `silent=false`, typically solvers will print information about the problem and their progress when solving it.

If you are experiencing performance issues, it is important to isolate which of these three steps is taking the most time. If everything is working as expected, step (1) should be the fastest, then step (2), then step (3). If solving the problem is much faster than formulating it, there is likely a performance problem.

## Avoid scalar indexing for faster problem formulation

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

