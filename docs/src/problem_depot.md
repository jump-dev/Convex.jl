Problem Depot
=============

Convex.jl has a submodule, `ProblemDepot` which holds a collection of convex optimization problems. The problems are used by Convex itself to test and benchmark its code, but can also be used by solvers to test and benchmark their code.

ProblemDepot has two main methods for accessing these problems: `Convex.ProblemDepot.run_tests` and `Convex.ProblemDepot.suite`.

For example, to test the solver SCS on all the problems of the depot except the mixed-integer problems (which it cannot handle), run
```julia
using Convex, SCS, Test
@testset "SCS" begin
    Convex.ProblemDepot.run_tests(; exclude=[r"mip"]) do p
        solve!(p, SCSSolver(verbose=0, eps=1e-6))
    end
end
```

## Reference

```@docs
Convex.ProblemDepot.run_tests
Convex.ProblemDepot.suite
Convex.ProblemDepot.PROBLEMS
```
