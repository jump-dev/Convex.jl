# Copyright (c) 2014: Madeleine Udell and contributors
#
# Use of this source code is governed by a BSD-style license that can be found
# in the LICENSE file or at https://opensource.org/license/bsd-2-clause

struct SparseAffineOperation{T}
    matrix::SPARSE_MATRIX{T}
    vector::SPARSE_VECTOR{T}
end

# function Base.convert(
#     ::Type{SparseAffineOperation{T}},
#     obj::SparseAffineOperation,
# ) where {T}
#     mat = convert(SparseMatrixCSC{T}, obj.matrix)
#     vec = convert(Vector{T}, obj.vector)
#     return SparseAffineOperation{T}(mat, vec)
# end

function SparseAffineOperation(
    A::AbstractMatrix{T},
    b::AbstractVector{T},
) where {T}
    return SparseAffineOperation{T}(create_sparse(T, A), SPARSE_VECTOR{T}(b))
end

struct LinkedSparseAffineOperation{T}
    previous::Union{LinkedSparseAffineOperation{T}, Nothing}
    current::SparseAffineOperation{T}
end

function LinkedSparseAffineOperation(sp::SparseAffineOperation{T}) where T
    return LinkedSparseAffineOperation{T}(nothing, sp)
end

function SparseAffineOperation(linked::LinkedSparseAffineOperation)
    if linked.previous === nothing
        return linked.current
    else
        return compose(linked.current, SparseAffineOperation(linked.previous))
    end
end

mutable struct SparseTape{T}
    operation::LinkedSparseAffineOperation{T}
    variables::Vector{MOI.VariableIndex}
    function SparseTape{T}(
        operation::SparseAffineOperation{T},
        variables::Vector{MOI.VariableIndex},
    ) where {T}
        return new(LinkedSparseAffineOperation(operation), variables)
    end

    function SparseTape(
        operation::SparseAffineOperation{T},
        variables::Vector{MOI.VariableIndex},
    ) where {T}
        return SparseTape{T}(operation, variables)
    end

    function SparseTape(
        operation::LinkedSparseAffineOperation{T},
        variables::Vector{MOI.VariableIndex},
    ) where {T}
        return new{T}(operation, variables)
    end
    function SparseTape{T}(
        operation::LinkedSparseAffineOperation{T},
        variables::Vector{MOI.VariableIndex},
    ) where {T}
        return new{T}(operation, variables)
    end
end

function Base.:(==)(tape1::SparseAffineOperation, tape2::SparseAffineOperation)
    return tape1.matrix == tape1.matrix && tape1.vector == tape2.vector
end
function Base.hash(tape::SparseAffineOperation, h::UInt)
    return hash(typeof(tape), hash(tape.matrix, hash(tape.vector, h)))
end

function Base.:(==)(tape1::SparseTape, tape2::SparseTape)
    return SparseAffineOperation(tape1) == SparseAffineOperation(tape2) &&
           tape1.variables == tape2.variables
end
function Base.hash(tape::SparseTape, h::UInt)
    return hash(typeof(tape), hash(SparseAffineOperation(tape), hash(tape.variables, h)))
end

MOI.output_dimension(v::SparseTape) = size(v.operation.current.matrix, 1)

function SparseAffineOperation(tape::SparseTape)# -> SparseAffineOperation
    return SparseAffineOperation(tape.operation)
end

function compose(A::SparseAffineOperation, B::SparseAffineOperation)
    vec = A.vector + A.matrix * B.vector
    mat = A.matrix * B.matrix
    return SparseAffineOperation(mat, vec)
end

#### SparseTape

function add_operation(tape::SparseTape{T}, op::SparseAffineOperation) where {T}
    return SparseTape(compose(op, SparseAffineOperation(tape)), tape.variables)
end

Base.real(tape::SparseTape) = tape

function Base.imag(c::SparseTape{T}) where {T}
    n = MOI.output_dimension(c)
    m = length(c.variables)
    mat = spzeros(T, n, m)
    v = spzeros(T, n)
    op = SparseAffineOperation(mat, v)

    # Hack re-use variables from input
    # I think this is OK bc the operation is all zeros
    return SparseTape(op, c.variables)
end
