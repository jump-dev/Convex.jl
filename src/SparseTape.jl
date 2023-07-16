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

# function SparseAffineOperation(
#     A::GBMatrix{T, T},
#     b::GBVector{T, T},
# ) where {T}
#     return SparseAffineOperation{T}(A, b)
# end

# function SparseAffineOperation(A::AbstractSparseMatrix, b)
#     return SparseAffineOperation(SparseMatrixCSC(A), b)
# end
# SparseAffineOperation(A, b) = SparseAffineOperation(create_sparse(A), b)

mutable struct SparseTape{T}
    operations::Vector{SparseAffineOperation{T}}
    variables::Vector{MOI.VariableIndex}
    function SparseTape{T}(
        operations::Vector{SparseAffineOperation{T}},
        variables::Vector{MOI.VariableIndex},
    ) where {T}
        # Is this necessary?
        # if !issorted(variables; by = x->x.value)
        #     p = sortperm(variables; by = x->x.value)
        #     op = foldl(compose, operations)
        #     matrix = op.matrix[:, p]
        #     vector = op.vector
        #     operations = [SparseAffineOperation(matrix, vector)]
        #     variables = variables[p]
        # end
        return new(operations, variables)
    end

    function SparseTape(
        operations::Vector{SparseAffineOperation{T}},
        variables::Vector{MOI.VariableIndex},
    ) where {T}
        return SparseTape{T}(operations, variables)
    end
end

MOI.output_dimension(v::SparseTape) = size(v.operations[1].matrix, 1)

function SparseAffineOperation(tape::SparseTape)# -> SparseAffineOperation
    return foldl(compose, tape.operations)
end

function compose(A::SparseAffineOperation, B::SparseAffineOperation)
    vec = A.vector + A.matrix * B.vector
    mat = A.matrix * B.matrix
    return SparseAffineOperation(mat, vec)
end

function collapse(sparse_tape::SparseTape) # -> SparseTape
    op = SparseAffineOperation(sparse_tape)
    return SparseTape([op], sparse_tape.variables)
end
#### SparseTape

function add_operation(tape::SparseTape{T}, op::SparseAffineOperation) where {T}
    tape2 = SparseTape(copy(tape.operations), tape.variables)
    pushfirst!(tape2.operations, op)
    return tape2::SparseTape{T}
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
    return SparseTape([op], c.variables)
end
