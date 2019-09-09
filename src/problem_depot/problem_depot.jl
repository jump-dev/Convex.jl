# Some code in `src/problem_depot` was modified from MathOptInterface
# which is available under an MIT license (see LICENSE).
module ProblemDepot
using BenchmarkTools, Test

using Convex

using Random
import LinearAlgebra.eigen
import LinearAlgebra.I
import LinearAlgebra.opnorm
import Random.shuffle
import Statistics.mean
using LinearAlgebra

const PROBLEMS = Dict{String, Dict{String, Function}}()

"""
    suite(
        handle_problem!::Function;
        exclude::Vector{Regex} = Regex[]
    )

Create a suite of benchmarks. `handle_problem!` should be a function that takes one
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
function suite(handle_problem!::Function, args...; exclude::Vector{Regex} = Regex[])
    group = BenchmarkGroup()
    for (class, dict) in PROBLEMS
        any(occursin.(exclude, Ref(class))) && continue
        for (name, func) in dict
            any(occursin.(exclude, Ref(name))) && continue
            group[name] = @benchmarkable $func($handle_problem!, args...)
        end
    end
    return group
end


function run_test(handle_problem!::Function; exclude::Vector{Regex} = Regex[], T=Float64, atol=1e-4, rtol=0.0, test = Val(true))
    for (class, dict) in PROBLEMS
        any(occursin.(exclude, Ref(class))) && continue
        @testset "$class" begin
            for (name, func) in dict
                any(occursin.(exclude, Ref(name))) && continue
                @testset "$name" begin
                    func(handle_problem!, test, atol, rtol, T)
                end
            end
        end
    end
end

###
### Benchmarks
###

macro add_problem(prefix, q)
    @assert prefix isa Symbol
    if q.head == :block
        f = q.args[2]
    elseif q.head == :function
        f = q
    else
        error("head $(q.head) unexpected")
    end
    name = f.args[1].args[1]
    if name isa Expr
        name = name.args[1]
    end
    return quote
        $(esc(f))
        dict = get!(PROBLEMS, String($(Base.Meta.quot(prefix))), Dict{String,Function}())
        dict[String($(Base.Meta.quot(name)))] = $(esc(name))
    end
end

eye(n, T) = Matrix{T}(I, n, n)
eye(n) = Matrix{Float64}(I, n, n)

include("problems/affine.jl")
include("problems/constant.jl")
include("problems/exp.jl")
include("problems/lp.jl")
include("problems/mip.jl")
include("problems/sdp_and_exp.jl")
include("problems/sdp.jl")
include("problems/socp.jl")

end
