# `BENCHMARKS`, `suite`, and `add_benchmark` were taken from MathOptInterface
# which is available under an MIT license (see LICENSE).
module Benchmarks

using BenchmarkTools
using Convex, LinearAlgebra

const BENCHMARKS = Dict{String, Function}()

"""
    suite(
        handle_problem::Function;
        exclude::Vector{Regex} = Regex[]
    )

Create a suite of benchmarks. `handle_problem` should be a function that takes one
argument, a Convex.jl `Problem` and processes it (e.g. `solve!` the problem with
a specific solver).

Use `exclude` to exclude a subset of benchmarks.

### Examples

```julia
suite() do p
    solve!(p, GLPK.Optimizer())
end
```
"""
function suite(handle_problem::Function; exclude::Vector{Regex} = Regex[])
    group = BenchmarkGroup()
    for (name, func) in BENCHMARKS
        any(occursin.(exclude, Ref(name))) && continue
        group[name] = @benchmarkable $func($handle_problem) setup=Convex.clearmemory()
    end
    return group
end

###
### Benchmarks
###

macro add_benchmark(f)
    name = f.args[1].args[1]
    return quote
        $(esc(f))
        BENCHMARKS[String($(Base.Meta.quot(name)))] = $(esc(name))
    end
end

eye(n) = Matrix(1.0I, n, n)
include("constraint_benchmarks.jl")
include("affine_benchmarks.jl")
include("sdp_benchmarks.jl")

end