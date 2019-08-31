# This file was copied from Transducers.jl
# which is available under an MIT license (see LICENSE).
using PkgBenchmark
include("pprinthelper.jl")
result = PkgBenchmark.readresults(joinpath(@__DIR__, "result.json"))
displayresult(result)
