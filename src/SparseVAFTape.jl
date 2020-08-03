# An alternate implementation of the VAF tape concept:
# Here, we represent each affine transformation by a single sparse matrix
# of uniform type. This helps avoid type issues, since every operation has the same type,
# and does not rely on correctly hitting optimized dispatches.
# On the other hand, we can't use optimized dispatches to help speed things up when
# we have extra structural information.

struct SparseAffineOperation{T <: AbstractSparseMatrix}
    matrix::T
end

SparseAffineOperation(A,b::AbstractSparseVector) = SparseAffineOperation(sparse(A), b)
SparseAffineOperation(A::AbstractSparseMatrix,b) = SparseAffineOperation(A, sparse(b))
SparseAffineOperation(A, b) = SparseAffineOperation(sparse(A), sparse(b))

function SparseAffineOperation(A::SparseMatrixCSC, b::SparseVector)
    T = eltype(A)
    # construct a sparse matrix representation of `Ax+b`
    # via `[A b; 0 1] * [x; 1] == [Ax+b; 1]`.
    # Thanks to Steve Diamond https://github.com/cvxgrp/cvxpy/issues/708#issuecomment-667692989
    n, m = size(A)
    A = sparse(A)
    mat = [A b; transpose(spzeros(T, m)) one(T)]
    SparseAffineOperation(mat)
end

struct SparseVAFTape{T}
    operations::Vector{SparseAffineOperation{T}}
    variables::Vector{MOI.VariableIndex}
end

MOI.output_dimension(v::SparseVAFTape) = size(v.operations[1].matrix, 1) - 1


_unwrap(a::SparseAffineOperation) = a.matrix
_unwrap(a::AbstractSparseArray) = a

function AffineOperation(sparse_tape::SparseVAFTape)
    if length(sparse_tape.operations) > 1
        mat = foldl((a,b) -> _unwrap(a) * _unwrap(b), sparse_tape.operations)
    else
        mat = only(sparse_tape.operations).matrix
    end
    b = mat[1:end-1, end]
    A = mat[1:end-1, 1:end-1]
    AffineOperation(A, b)
end


#### SparseVAFTape

function add_operation!(tape::SparseVAFTape, op::SparseAffineOperation)
    pushfirst!(tape.operations, op)
    tape
end

# the lack of `!` isn't quite right, since we mutate the `tape`,
# but we shouldn't run into problems since outside packages aren't using `VAFTape`'s.
function MOIU.operate(::typeof(-), ::Type{T},
    tape::SparseVAFTape) where {T}
    d = MOI.output_dimension(tape)
    return add_operation!(tape, SparseAffineOperation(sparse(-one(T)*I, d, d), Zero(d)))
end

function MOIU.operate(::typeof(+), ::Type{T}, v::AbstractVector,
                            tape::SparseVAFTape) where {T}
    d = length(v)
    return add_operation!(tape, SparseAffineOperation(sparse(one(T)*I, d, d), v))
end

function MOIU.operate(::typeof(-), ::Type{T}, tape::SparseVAFTape,
    v::AbstractVector) where {T}
    d = length(v)
    return add_operation!(tape, SparseAffineOperation(sparse(one(T)*I, d, d), -v))
end

function MOIU.operate(::typeof(-), ::Type{T},
                    v::AbstractVector, tape::SparseVAFTape) where {T}
    d = length(v)
    return add_operation!(tape, SparseAffineOperation(sparse(-one(T)*I, d, d), v))
end

function MOIU.operate(::typeof(*), ::Type{T}, A::AbstractMatrix,
    tape::SparseVAFTape) where {T}
    return add_operation!(tape, SparseAffineOperation(A, Zero(size(A,1))))
end

function MOIU.operate(::typeof(sum), ::Type{T}, tape::SparseVAFTape) where {T}
    d = MOI.output_dimension(tape)
    # doesn't seem ideal for a sparse representation...
    A = ones(T, 1, d) 
    return add_operation!(tape, SparseAffineOperation(A, Zero(size(A,1))))
end
