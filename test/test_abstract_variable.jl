# We provide a non-`Variable` implementation of the `AbstractVariable` interface
# to test that only the interface is used (and not, e.g. direct field access).
    
module TypedVectors
using Convex

# To make sure `Convex` isn't using field access on `AbstractVariable`'s
# we'll use a global dictionary to store information about each instance
# our of mock variable type, `TypedVector`.
const global_cache = Dict{UInt64, Any}()

mutable struct TypedVector{T} <: Convex.AbstractVariable
    head::Symbol
    id_hash::UInt64
    size::Tuple{Int, Int}
    function TypedVector{T}(d) where {T}
        this = new(:ConstSizeVariable, 0, (d,1))
        this.id_hash = objectid(this)
        Convex.id_to_variables[this.id_hash] = this
        global_cache[this.id_hash] = Dict(  :value => nothing,
                                            :sign => T <: Complex ? ComplexSign() : NoSign(), 
                                            :vartype => ContVar,
                                            :constraints => Constraint[],
                                            :vexity => AffineVexity())
        this
    end
end

Convex.value(x::TypedVector) = global_cache[x.id_hash][:value]

Convex.value!(x::TypedVector, v::AbstractArray) = global_cache[x.id_hash][:value] = v
Convex.value!(x::TypedVector, v::Number) = global_cache[x.id_hash][:value] = v

Convex.vexity(x::TypedVector) = global_cache[x.id_hash][:vexity]
Convex.vexity!(x::TypedVector, v::Vexity) = global_cache[x.id_hash][:vexity] = v

Convex.sign(x::TypedVector) = global_cache[x.id_hash][:sign]
Convex.sign!(x::TypedVector, s::Sign) = global_cache[x.id_hash][:sign] = s

Convex.vartype(x::TypedVector) = global_cache[x.id_hash][:vartype]
Convex.vartype!(x::TypedVector, s::Convex.VarType) = global_cache[x.id_hash][:vartype] = s

Convex.constraints(x::TypedVector) = global_cache[x.id_hash][:constraints]
Convex.add_constraint!(x::TypedVector, s::Constraint) = push!(global_cache[x.id_hash][:constraints], s)

Convex.eltype(x::TypedVector{T}) where {T} = T

end

import .TypedVectors

@testset "AbstractVariable interface: $solver" for solver in solvers
    # Let us solve a basic problem from `test_affine.jl`

    x = TypedVectors.TypedVector{BigFloat}(1)
    y = TypedVectors.TypedVector{BigFloat}(1)
    p = minimize(x + y, [x >= 3, y >= 2])
    @test vexity(p) == AffineVexity()
    solve!(p, solver)
    @test p.optval ≈ 5 atol=TOL
    @test evaluate(x + y) ≈ 5 atol=TOL
    @test Convex.eltype(x) == BigFloat

    add_constraint!(x, x >= 4)
    solve!(p, solver)
    @test p.optval ≈ 6 atol=TOL
    @test evaluate(x + y) ≈ 6 atol=TOL
    @test length(constraints(x)) == 1

end


# Let us do another example of custom variable types, but using field access for simplicity
  
module DensityMatricies
using Convex

mutable struct DensityMatrix <: Convex.AbstractVariable
    head::Symbol
    id_hash::UInt64
    size::Tuple{Int, Int}
    value::Convex.ValueOrNothing
    vexity::Vexity
    function DensityMatrix(d)
        this = new(:DensityMatrix, 0, (d,d), nothing, Convex.AffineVexity())
        this.id_hash = objectid(this)
        Convex.id_to_variables[this.id_hash] = this
        this
    end
end
Convex.constraints(ρ::DensityMatrix) = [ ρ ⪰ 0, tr(ρ) == 1 ]
Convex.sign(::DensityMatrix) = Convex.ComplexSign()
Convex.vartype(::DensityMatrix) = Convex.ContVar
Convex.eltype(::DensityMatrix) = ComplexF64

end

import .DensityMatricies
import LinearAlgebra
@testset "DensityMatrix custom variable type: $solver" for solver in solvers

    if can_solve_sdp(solver)
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
        solve!(prob, solver)

        @test prob.optval ≈ e_val atol = TOL
        @test evaluate(ρ) ≈ proj atol = TOL
    end
end
