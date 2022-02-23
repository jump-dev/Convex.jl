# We provide a non-`Variable` implementation of the `AbstractVariable` interface
# to test that only the interface is used (and not, e.g. direct field access).
solver = MOI.OptimizerWithAttributes(SCS.Optimizer, "verbose" => 0)
TOL = 1e-3
module DictVectors
using Convex

# To make sure `Convex` isn't using field access on `AbstractVariable`'s
# we'll use a global dictionary to store information about each instance
# our of mock variable type, `DictVector`.
const global_cache = Dict{UInt64, Any}()

mutable struct DictVector{T} <: Convex.AbstractVariable
    head::Symbol
    id_hash::UInt64
    size::Tuple{Int, Int}
    function DictVector{T}(d) where {T}
        this = new(:DictVector, 0, (d,1))
        this.id_hash = objectid(this)
        global_cache[this.id_hash] = Dict(  :value => nothing,
                                            :sign => T <: Complex ? ComplexSign() : NoSign(),
                                            :vartype => ContVar,
                                            :constraints => Constraint[],
                                            :vexity => Convex.AffineVexity())
        this
    end
end

Convex.evaluate(x::DictVector) = global_cache[x.id_hash][:value]

Convex.set_value!(x::DictVector, v::AbstractArray) = global_cache[x.id_hash][:value] = v
Convex.set_value!(x::DictVector, v::Number) = global_cache[x.id_hash][:value] = v

Convex.vexity(x::DictVector) = global_cache[x.id_hash][:vexity]
Convex.vexity!(x::DictVector, v::Convex.Vexity) = global_cache[x.id_hash][:vexity] = v

Convex.sign(x::DictVector) = global_cache[x.id_hash][:sign]
Convex.sign!(x::DictVector, s::Convex.Sign) = global_cache[x.id_hash][:sign] = s

Convex.vartype(x::DictVector) = global_cache[x.id_hash][:vartype]
Convex.vartype!(x::DictVector, s::Convex.VarType) = global_cache[x.id_hash][:vartype] = s

Convex.constraints(x::DictVector) = global_cache[x.id_hash][:constraints]
Convex.add_constraint!(x::DictVector, s::Convex.Constraint) = push!(global_cache[x.id_hash][:constraints], s)

end

import .DictVectors

@testset "AbstractVariable interface" begin
    # Let us solve a basic problem from `test_affine.jl`

    x = DictVectors.DictVector{BigFloat}(1)
    y = DictVectors.DictVector{BigFloat}(1)
    p = minimize(x + y, [x >= 3, y >= 2])
    @test vexity(p) == Convex.AffineVexity()
    solve!(p, solver())
    @test p.optval ≈ 5 atol=TOL
    @test evaluate(x + y) ≈ 5 atol=TOL

    add_constraint!(x, x >= 4)
    solve!(p, solver())
    @test p.optval ≈ 6 atol=TOL
    @test evaluate(x + y) ≈ 6 atol=TOL
    @test length(constraints(x)) == 1

end


# Let us do another example of custom variable types, but using field access for simplicity

module DensityMatricies
using Convex

mutable struct DensityMatrix{T} <: Convex.AbstractVariable
    head::Symbol
    id_hash::UInt64
    size::Tuple{Int, Int}
    value::Convex.ValueOrNothing
    vexity::Convex.Vexity
    function DensityMatrix(d)
        this = new{ComplexF64}(:DensityMatrix, 0, (d,d), nothing, Convex.AffineVexity())
        this.id_hash = objectid(this)
        this
    end
end
Convex.constraints(ρ::DensityMatrix) = [ ρ ⪰ 0, tr(ρ) == 1 ]
Convex.sign(::DensityMatrix) = Convex.ComplexSign()
Convex.vartype(::DensityMatrix) = Convex.ContVar

end

import .DensityMatricies
import LinearAlgebra
@testset "DensityMatrix custom variable type" begin
    X = randn(4,4) + im*rand(4,4); X = X+X'
    # `X` is Hermitian and non-degenerate (with probability 1)
    # Let us calculate the projection onto the eigenspace associated
    # to the maximum eigenvalue
    e_vals, e_vecs = LinearAlgebra.eigen(LinearAlgebra.Hermitian(X))
    e_val, idx = findmax(e_vals)
    e_vec = e_vecs[:, idx]
    proj = e_vec * e_vec'

    # found it! Now let us do it again via an SDP
    ρ = DensityMatricies.DensityMatrix(4)

    prob = maximize( real(tr(ρ*X)) )
    solve!(prob, solver())

    @test prob.optval ≈ e_val atol = TOL
    @test evaluate(ρ) ≈ proj atol = TOL
end


module ProbabilityVectors
using Convex
mutable struct ProbabilityVector <: Convex.AbstractVariable
    head::Symbol
    id_hash::UInt64
    size::Tuple{Int, Int}
    value::Convex.ValueOrNothing
    vexity::Convex.Vexity
    function ProbabilityVector(d)
        this = new(:ProbabilityVector, 0, (d,1), nothing, Convex.AffineVexity())
        this.id_hash = objectid(this)
        this
    end
end
Convex.constraints(p::ProbabilityVector) = [ sum(p) == 1 ]
Convex.sign(::ProbabilityVector) = Convex.Positive()
Convex.vartype(::ProbabilityVector) = Convex.ContVar

(p::ProbabilityVector)(x) = dot(p, x)

end

using .ProbabilityVectors

@testset "ProbabilityVectors" begin
    p = ProbabilityVectors.ProbabilityVector(3)
    x = [1.0, 2.0, 3.0]

    @test p(x) isa AbstractExpr
    @test sign(p) == Positive()
    prob = minimize( p(x) )
    solve!(prob, solver())
    @test prob.optval ≈ 1.0 atol=TOL
    @test evaluate(p(x)) ≈ 1.0 atol=TOL
    @test evaluate(p) ≈ [1.0, 0.0, 0.0] atol=TOL

end
