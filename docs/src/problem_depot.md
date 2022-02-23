Problem Depot
=============

Convex.jl has a submodule, `ProblemDepot` which holds a collection of convex optimization problems. The problems are used by Convex itself to test and benchmark its code, but can also be used by solvers to test and benchmark their code. These tests have been used with many solvers at [ConvexTests.jl](https://github.com/ericphanson/ConvexTests.jl).

ProblemDepot has two main methods for accessing these problems: `Convex.ProblemDepot.run_tests` and `Convex.ProblemDepot.benchmark_suite`.

For example, to test the solver SCS on all the problems of the depot except the mixed-integer problems (which it cannot handle), run

```julia
using Convex, SCS, Test
@testset "SCS" begin
    Convex.ProblemDepot.run_tests(; exclude=[r"mip"]) do p
        solve!(p, MOI.OptimizerWithAttributes(SCS.Optimizer, "verbose" => 0, "eps_abs" => 1e-6))
    end
end
```

How to write a ProblemDepot problem
-----------------------------------

The problems are organized into folders in `src/problem_depot/problems`. Each is written as a function, annotated by `@add_problem`, and a name, which is used to group the problems. For example, here is a simple problem:

```julia
@add_problem affine function affine_negate_atom(handle_problem!, ::Val{test}, atol, rtol, ::Type{T}) where {T, test}
    x = Variable()
    p = minimize(-x, [x <= 0])
    if test
        @test vexity(p) == AffineVexity()
    end
    handle_problem!(p)
    if test
        @test p.optval ≈ 0 atol=atol rtol=rtol
        @test evaluate(-x) ≈ 0 atol=atol rtol=rtol
    end
end
```

The `@add_problem` call adds the problem to the registry of problems in [`Convex.ProblemDepot.PROBLEMS`](@ref), which in turn is used by [`Convex.ProblemDepot.run_tests`](@ref) and [`Convex.ProblemDepot.benchmark_suite`](@ref). Next, `affine` is the grouping of the problem; this problem came from one of the affine tests, and in particular is testing the negation atom. Next is the function signature:

```julia
function affine_negate_atom(handle_problem!, ::Val{test}, atol, rtol, ::Type{T}) where {T, test}
```

this should be the same for every problem, except for the name, which is a description of the problem. It should include what kind of atoms it uses (`affine` in this case), so that certain kinds of atoms can be ruled out by the `exclude` keyword to [`Convex.ProblemDepot.run_tests`](@ref) and [`Convex.ProblemDepot.benchmark_suite`](@ref); for example, many solvers cannot solve mixed-integer problems, so `mip` is included in the name of such problems.

Then begins the body of the problem. It is setup like any other Convex.jl problem, only `handle_problem!` is called instead of `solve!`. This allows particular solvers to be used (via e.g. choosing `handle_problem! = p -> solve!(p, solver)`), or for any other function of the problem. Tests should be included and gated behind `if test` blocks, so that tests can be skipped for benchmarking, or in the case that the problem is not in fact solved during `handle_problem!`.

The fact that the problems may not be solved during `handle_problem!` brings with it a small complication: any command that assumes the problem has been solved should be behind an `if test` check. For example, in some of the problems, `real(evaluate(x))` is used, for a variable `x`; perhaps as

```julia
x_re = real(evaluate(x))
if test
    @test x_re = ...
end
```

However, if the problem `x` is used in has not been solved, then `evaluate(x) === nothing`, and `real(nothing)` throws an error. So instead, this should be rewritten as

```julia
if test
    x_re = real(evaluate(x))
    @test x_re = ...
end
```

Benchmark-only problems
-----------------------

To add problems for benchmarking without tests, place problems in `src/problem_depot/problems/benchmark`, and include `benchmark` in the name. These problems will be automatically skipped during `run_tests` calls. For example, to benchmark the time it takes to add an SDP constraint, we have the problem

```julia
@add_problem constraints_benchmark function sdp_constraint(handle_problem!, args...)
    p = satisfy()
    x = Variable(44, 44) # 990 vectorized entries
    push!(p.constraints, x ⪰ 0)
    handle_problem!(p)
    nothing
end
```

However, this "problem" has no tests or interesting content for testing, so we skip it during testing.
Note, we use `args...` in the function signature so that it may be called with the standard function signature

```julia
f(handle_problem!, ::Val{test}, atol, rtol, ::Type{T}) where {T, test}
```

Reference
---------

```@docs
Convex.ProblemDepot.run_tests
Convex.ProblemDepot.benchmark_suite
Convex.ProblemDepot.foreach_problem
Convex.ProblemDepot.PROBLEMS
```
