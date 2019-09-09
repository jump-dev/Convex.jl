using Pkg
tempdir = mktempdir()
Pkg.activate(tempdir)
Pkg.develop(PackageSpec(path=joinpath(@__DIR__, "..")))
Pkg.add(["BenchmarkTools", "SCS", "ECOS", "PkgBenchmark"])
Pkg.resolve()

using Convex: Convex, ProblemDepot
using SCS: SCSSolver
using ECOS: ECOSSolver
using BenchmarkTools


const SUITE = BenchmarkGroup()

# SUITE["SCS"] = ProblemDepot.benchmark_suite(p -> Convex.solve!(p, SCSSolver(verbose=0)))
# SUITE["ECOS"] = ProblemDepot.benchmark_suite(p -> Convex.solve!(p, ECOSSolver(verbose=0));
    # exclude = [r"sdp", r"exp"])
SUITE["formulation"] = ProblemDepot.benchmark_suite(Convex.conic_problem; exclude=[r"fix"])
