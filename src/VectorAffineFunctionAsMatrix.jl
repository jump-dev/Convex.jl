# # A simple type representing a vector of zeros. Maybe should include the size or use FillArrays or similar.
# struct Zero
#     len::Int
# end

# Base.:(+)(a, ::Zero) = a
# Base.:(+)(::Zero, a) = a
# Base.:(+)(z1::Zero, z2::Zero) = (@assert z1.len == z2.len; z1)
# Base.:(-)(::Zero, a) = -a
# Base.:(-)(a, ::Zero) = a
# Base.:(-)(z::Zero) = z
# Base.:(*)(A, z::Zero) = Zero(size(A, 1))
# Base.:(*)(::LinearAlgebra.UniformScaling, z::Zero) = z
# Base.size(z::Zero) = (z.len,)
# Base.length(z::Zero) = z.len
# SparseArrays.sparse(z::Zero) = spzeros(z.len)
# Base.convert(::Type{Vector{T}}, z::Zero) where {T} = zeros(T, z.len)

# This is a variant of MOI.VectorAffineFunction which represents
# the transformation `matrix * variables + vector` lazily.
struct VectorAffineFunctionAsMatrix{T,V}
    aff::SparseAffineOperation{T}
    variables::V
end

function Base.isequal(
    A::VectorAffineFunctionAsMatrix,
    B::VectorAffineFunctionAsMatrix,
)
    return isequal(A.variables, B.variables) && isequal(A.aff, B.aff)
end

function MOI.output_dimension(v::VectorAffineFunctionAsMatrix)
    return size(v.matrix, 1)
end

function to_vaf(tape::SparseTape)
    op = SparseAffineOperation(tape)
    return to_vaf(VectorAffineFunctionAsMatrix(op, tape.variables))
end

# convert to a usual VAF
function to_vaf(vaf_as_matrix::VectorAffineFunctionAsMatrix{T}) where {T}
    I, J, V = findnz(vaf_as_matrix.aff.matrix)
    vats = MOI.VectorAffineTerm{T}[]
    for n in eachindex(I, J, V)
        i = I[n]
        j = J[n]
        v = V[n]
        push!(
            vats,
            MOI.VectorAffineTerm{T}(
                i,
                MOI.ScalarAffineTerm{T}(v, vaf_as_matrix.variables[j]),
            ),
        )
    end

    return MOI.VectorAffineFunction{T}(vats, vaf_as_matrix.aff.vector)
end

# method for adding constraints and coverting to standard VAFs as needed
function MOI_add_constraint(model, f, set)
    return MOI.add_constraint(model, f, set)
end

function MOI_add_constraint(
    model,
    f::Union{VectorAffineFunctionAsMatrix,SparseTape},
    set,
)
    return MOI_add_constraint(model, to_vaf(f), set)
end
