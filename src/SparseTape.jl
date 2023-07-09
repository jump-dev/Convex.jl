struct SparseAffineOperation{T}
    matrix::SparseMatrixCSC{T}
    vector::Vector{T}
end

function Base.convert(
    ::Type{SparseAffineOperation{T}},
    obj::SparseAffineOperation,
) where {T}
    mat = convert(SparseMatrixCSC{T}, obj.matrix)
    vec = convert(Vector{T}, obj.vector)
    return SparseAffineOperation{T}(mat, vec)
end

function SparseAffineOperation(
    A::SparseMatrixCSC{T},
    b::AbstractVector,
) where {T}
    return SparseAffineOperation(A, collect(T, b))
end

function SparseAffineOperation(A::AbstractSparseMatrix, b)
    return SparseAffineOperation(SparseMatrixCSC(A), b)
end
SparseAffineOperation(A, b) = SparseAffineOperation(sparse(A), b)

struct SparseTape{T}
    operations::Vector{SparseAffineOperation{T}}
    variables::Vector{MOI.VariableIndex}
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
    v = zeros(T, n)
    op = SparseAffineOperation(mat, v)

    # Hack re-use variables from input
    # I think this is OK bc the operation is all zeros
    return SparseTape([op], c.variables)
end
