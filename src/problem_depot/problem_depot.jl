# Some code in `src/problem_depot` was modified from MathOptInterface
# which is available under an MIT license (see LICENSE).

module ProblemDepot
using BenchmarkTools, Test
using MathOptInterface
const MOI = MathOptInterface
using Convex
using LinearAlgebra
using LinearAlgebra: eigen, I, opnorm
using Convex: AffineVexity, ConcaveVexity, ConvexVexity

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

See also [`run_tests`](@ref) and [`benchmark_suite`](@ref) for helpers
to use these problems in testing or benchmarking.

### Examples

```julia
julia> PROBLEMS["affine"]["affine_diag_atom"]
affine_diag_atom (generic function with 1 method)
```
"""
const PROBLEMS = Dict{String, Dict{String, Function}}()

"""
    foreach_problem(apply::Function, [class::String],
        problems::Union{Nothing, Vector{String}, Vector{Regex}} = nothing;
        exclude::Vector{Regex} = Regex[])

Provides a convience method for iterating over problems in [`PROBLEMS`](@ref).
For each problem in [`PROBLEMS`](@ref), apply the function `apply`, which
takes two arguments: the name of the function associated to the problem,
and the function associated to the problem itself.

Optionally, pass a second argument `class` to only iterate over a class of
problems (`class` should satsify `class ∈ keys(PROBLEMS)`), and pass third
argument `problems` to only allow certain problems (specified by exact names or
regex). Use the `exclude` keyword argument to exclude problems by regex.
"""
function foreach_problem(   apply::Function,
                            problems::Union{Nothing, Vector{String}, Vector{Regex}} = nothing;
                            exclude::Vector{Regex} = Regex[])
    for class in keys(PROBLEMS)
        any(occursin.(exclude, Ref(class))) && continue
        foreach_problem(apply, class, problems; exclude=exclude)
    end
end

function foreach_problem(   apply::Function,
                            class::String,
                            problems::Union{Nothing, Vector{String}, Vector{Regex}} = nothing;
                            exclude::Vector{Regex} = Regex[])
    for (name, func) in PROBLEMS[class]
        any(occursin.(exclude, Ref(name))) && continue
        if problems !== nothing
            problems isa Vector{String} && !(name ∈ problems) && continue
            problems isa Vector{Regex} && !any(occursin.(problems, Ref(name))) && continue
        end
        apply(name, func)
    end
end


"""
    run_tests(
        handle_problem!::Function;
        problems::Union{Nothing, Vector{String}, Vector{Regex}} = nothing;
        exclude::Vector{Regex} = Regex[],
        T=Float64, atol=1e-3, rtol=0.0,
    )

Run a set of tests. `handle_problem!` should be a function that takes one
argument, a Convex.jl `Problem` and processes it (e.g. `solve!` the problem with
a specific solver).

Use `exclude` to exclude a subset of sets; automatically excludes
`r"benchmark"`. Optionally, pass a second argument `problems` to only allow certain problems
(specified by exact names or regex). The test tolerances specified by `atol` and
`rtol`. Set `T` to choose a numeric type for the problem. Currently
this is only used for choosing the type parameter of the underlying
MathOptInterface model, but not for the actual problem data.

### Examples

```julia
run_tests(exclude=[r"mip"]) do p
    solve!(p, () -> SCS.Optimizer(verbose=0))
end
```
"""
function run_tests( handle_problem!::Function,
                    problems::Union{Nothing, Vector{String}, Vector{Regex}} = nothing;
                    exclude::Vector{Regex} = Regex[], T=Float64, atol=1e-3, rtol=0.0)
    push!(exclude, r"benchmark")
    for class in keys(PROBLEMS)
        any(occursin.(exclude, Ref(class))) && continue
        @testset "$class" begin
            foreach_problem(class, problems; exclude=exclude) do name, problem_func
                @testset "$name" begin
                    problem_func(handle_problem!, Val(true), atol, rtol, T)
                end
            end
        end
    end
end


"""
    benchmark_suite(
        handle_problem!::Function,
        problems::Union{Nothing, Vector{String}, Vector{Regex}} = nothing;
        exclude::Vector{Regex} = Regex[],
        test = Val(false),
        T=Float64, atol=1e-3, rtol=0.0,
    )

Create a benchmark_suite of benchmarks. `handle_problem!` should be a function
that takes one argument, a Convex.jl `Problem` and processes it (e.g. `solve!`
the problem with a specific solver). Pass a second argument `problems` to specify
run benchmarks only with certain problems (specified by exact names or regex).

Use `exclude` to exclude a subset of benchmarks. Optionally, pass a second
argument `problems` to only allow certain problems (specified by exact names or
regex). Set `test=true` to also check the answers, with tolerances specified by
`atol` and `rtol`. Set `T` to choose a numeric type for the problem. Currently
this is only used for choosing the type parameter of the underlying
MathOptInterface model, but not for the actual problem data.

### Examples

```julia
benchmark_suite(exclude=[r"mip"]) do p
    solve!(p, () -> SCS.Optimizer(verbose=0))
end
```
"""
function benchmark_suite(handle_problem!::Function,
                         problems::Union{Nothing, Vector{String}, Vector{Regex}} = nothing;
                         exclude::Vector{Regex} = Regex[],
                         T=Float64, atol=1e-3, rtol=0.0, test = Val(false))
    group = BenchmarkGroup()
    for class in keys(PROBLEMS)
        any(occursin.(exclude, Ref(class))) && continue
        group[class] = BenchmarkGroup()
        foreach_problem(class, problems; exclude=exclude) do name, problem_func
            group[class][name] = @benchmarkable $problem_func($handle_problem!, $test, $atol, $rtol, $T)
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
