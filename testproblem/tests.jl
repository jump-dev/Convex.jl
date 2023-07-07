using Convex, SCS, Test, Clarabel
using Convex.ProblemDepot: run_tests
import MathOptInterface as MOI
Convex.USE_SPARSE2() = true

@testset "SCS without duals" begin
    run_tests(
        exclude = [
            r"mip",
            # TODO(odow): investigate
            r"sdp_lieb_ando",
            # Tolerance issue with SCS 3.0
            r"sdp_Real_Variables_with_complex_equality_constraints",
            r"dual",
            r"constant",
        ],
    ) do p
        return solve!(
            p,
            MOI.OptimizerWithAttributes(
                SCS.Optimizer,
                "verbose" => 0,
                "eps_rel" => 1e-6,
                "eps_abs" => 1e-6,
            );
        )
    end
end

@testset "SCS" begin
    run_tests([r"dual"], exclude = [r"socp"]) do p
        return solve2!(
            p,
            MOI.OptimizerWithAttributes(
                SCS.Optimizer,
                "verbose" => 0,
                "eps_rel" => 1e-6,
                "eps_abs" => 1e-6,
            );
        )
    end
end

@testset "Clarabel" begin
    run_tests([r"dual"], exclude = [r"socp"]) do p
        return solve2!(p, Clarabel.Optimizer; silent_solver = true)
    end
end
