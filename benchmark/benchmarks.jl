using Pkg
tempdir = mktempdir()
Pkg.activate(tempdir)
Pkg.develop(PackageSpec(path=joinpath(@__DIR__, "..")))
Pkg.add(["BenchmarkTools", "SCS", "ECOS", "PkgBenchmark"])
Pkg.resolve()

using Convex: Convex
using SCS: SCSSolver
using ECOS: ECOSSolver
using BenchmarkTools


include("benchmarks/benchmarks.jl") # defines module Benchmarks

const SUITE = BenchmarkGroup()

SUITE["SCS"] = Benchmarks.suite(p -> Convex.solve!(p, SCSSolver(verbose=0)))
SUITE["ECOS"] = Benchmarks.suite(p -> Convex.solve!(p, ECOSSolver(verbose=0));
    exclude = [r"sdp", r"SDP"])
SUITE["formulation"] = Benchmarks.suite(Convex.conic_problem)
