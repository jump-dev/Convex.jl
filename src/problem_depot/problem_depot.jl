# Some code in `src/problem_depot` was modified from MathOptInterface
# which is available under an MIT license (see LICENSE).
module ProblemDepot
using BenchmarkTools, Test

using Convex
using LinearAlgebra
using LinearAlgebra: eigen, I, opnorm

randperm(d) = sortperm(rand(d))
shuffle(x) = x[randperm(length(x))]
mean(x) = sum(x) / length(x)
eye(n, T) = Matrix{T}(I, n, n)
eye(n) = Matrix{Float64}(I, n, n)

"""
    const PROBLEMS = Dict{String, Dict{String, Function}}()

A "depot" of Convex.jl problems, subdivided into categories.
Each problem is stored as a function with the signature
    
    f(handle_problem!, ::Val{test}, atol, rtol, ::Type{T}) where {T, test}
    
where `handle_problem!` specifies what to do with the `Problem` instance
(e.g., `solve!` it with a chosen solver), an option `test` to choose
whether or not to test the values (assuming it has been solved),
tolerances for the tests, and a numeric type in which the problem
should be specified (currently, this is not respected and all
problems are specified in `Float64` precision).

See also [`run_tests`](@ref) and [`suite`](@ref) for helpers
to use these problems in testing or benchmarking.

### Examples

```julia
julia> PROBLEMS["affine"]["affine_diag_atom"]
affine_diag_atom (generic function with 1 method)

```
"""
const PROBLEMS = Dict{String, Dict{String, Function}}()


"""
    run_tests(
        handle_problem!::Function;
        exclude::Vector{Regex} = Regex[],
        T=Float64, atol=1e-3, rtol=0.0, 
    )

Run a set of tests benchmarks. `handle_problem!` should be a function that takes one
argument, a Convex.jl `Problem` and processes it (e.g. `solve!` the problem with
a specific solver).

Use `exclude` to exclude a subset of sets. The test tolerances specified by `atol` and `rtol`.
Set `T` to choose a numeric type for the problem. Currently this option is not respected
and all problems are specified Float64` precision.

### Examples

```julia
run_tests(exclude=[r"mip"]) do p
    solve!(p, SCSSolver(verbose=0))
end
```
"""
function run_tests(handle_problem!::Function; exclude::Vector{Regex} = Regex[], T=Float64, atol=1e-3, rtol=0.0)
    for (class, dict) in PROBLEMS
        any(occursin.(exclude, Ref(class))) && continue
        @testset "$class" begin
            for (name, func) in dict
                any(occursin.(exclude, Ref(name))) && continue
                @testset "$name" begin
                    func(handle_problem!, Val(true), atol, rtol, T)
                end
            end
        end
    end
end


"""
    suite(
        handle_problem!::Function;
        exclude::Vector{Regex} = Regex[],
        test = Val(false),
        T=Float64, atol=1e-3, rtol=0.0, 
    )

Create a suite of benchmarks. `handle_problem!` should be a function that takes one
argument, a Convex.jl `Problem` and processes it (e.g. `solve!` the problem with
a specific solver).

Use `exclude` to exclude a subset of benchmarks. Set `test=true` to also check the
answers, with tolerances specified by `atol` and `rtol`. Set `T` to choose a numeric
type for the problem. Currently this option is not respected and all problems are
specified Float64` precision.

### Examples

```julia
suite(exclude=[r"mip"]) do p
    solve!(p, SCSSolver(verbose=0))
end
```
"""
function suite(handle_problem!::Function; exclude::Vector{Regex} = Regex[], T=Float64, atol=1e-3, rtol=0.0, test = Val(false))
    group = BenchmarkGroup()
    for (class, dict) in ProblemDepot.PROBLEMS
        any(occursin.(exclude, Ref(class))) && continue
        group[class] = class_group = BenchmarkGroup()
        for (name, func) in dict
            any(occursin.(exclude, Ref(name))) && continue
                class_group[name] = @benchmarkable $func($handle_problem!, $test, $atol, $rtol, $T)
        end
    end
    return group
end


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

include("problems/affine.jl")
include("problems/constant.jl")
include("problems/exp.jl")
include("problems/lp.jl")
include("problems/mip.jl")
include("problems/sdp_and_exp.jl")
include("problems/sdp.jl")
include("problems/socp.jl")

end
