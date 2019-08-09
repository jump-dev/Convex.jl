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
Convex.constraints!(x::TypedVector, s::Vector{Constraint}) = global_cache[x.id_hash][:constraints] = s

Convex.eltype(x::TypedVector{T}) where {T} = T

end

import .TypedVectors

@testset "AbstractVariable interface" for solver in solvers
    # Basic problem

    x = TypedVectors.TypedVector{BigFloat}(1)
    y = TypedVectors.TypedVector{BigFloat}(1)
    p = minimize(x + y, [x >= 3, y >= 2])
    @test vexity(p) == AffineVexity()
    solve!(p, solver)
    @test p.optval ≈ 5 atol=TOL
    @test evaluate(x + y) ≈ 5 atol=TOL
    @test Convex.eltype(x) == BigFloat
end