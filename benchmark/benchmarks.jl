using Pkg
tempdir = mktempdir()
Pkg.activate(tempdir)
Pkg.develop(PackageSpec(path = joinpath(@__DIR__, "..")))
Pkg.add(["BenchmarkTools", "PkgBenchmark", "MathOptInterface"])
Pkg.resolve()

using Convex: Convex, ProblemDepot
using BenchmarkTools
using MathOptInterface
const MOI = MathOptInterface
const MOIU = MOI.Utilities

const SUITE = BenchmarkGroup()

problems = [
    "constant_fix!_with_complex_numbers",
    "affine_dot_multiply_atom",
    "affine_hcat_atom",
    "affine_trace_atom",
    "exp_entropy_atom",
    "exp_log_perspective_atom",
    "socp_norm_2_atom",
    "socp_quad_form_atom",
    "socp_sum_squares_atom",
    "lp_norm_inf_atom",
    "lp_maximum_atom",
    "sdp_and_exp_log_det_atom",
    "sdp_norm2_atom",
    "sdp_lambda_min_atom",
    "sdp_sum_largest_eigs",
    "mip_integer_variables",
]

SUITE["formulation"] = ProblemDepot.benchmark_suite(problems) do problem
    opt = MOIU.MockOptimizer(MOIU.Model{Float64}())
    return Convex.Context(problem, opt)
end
