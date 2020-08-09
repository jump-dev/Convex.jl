using Convex, SCS, Test
using Convex.ProblemDepot: run_tests

Convex.USE_SPARSE() = false
Convex.USE_SPARSE() = true
Convex.USE_SPARSE3() = false

@testset "SCS" begin
    run_tests([r"dual"], exclude=[r"socp"]) do p
        solve2!(p, () -> SCS.Optimizer(verbose=0, eps=1e-6))
    end
end
