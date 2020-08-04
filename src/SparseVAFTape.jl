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

function MOIU.operate(::typeof(+), ::Type{T}, tape::SparseVAFTape, vs::AbstractVector...) where {T}
    v = sum(vs)
    d = length(v)
    return add_operation!(tape, SparseAffineOperation(sparse(one(T)*I, d, d), v))
end

function MOIU.operate(::typeof(+), ::Type{T}, v::AbstractVector, tape::SparseVAFTape, vs::AbstractVector...) where {T}
    return MOIU.operate(+, T, tape, v + sum(vs))
end


function MOIU.operate(::typeof(+), ::Type{T}, v::AbstractVector, tape::SparseVAFTape) where {T}
    return MOIU.operate(+, T, tape, v)
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


function MOIU.operate(::typeof(+), ::Type{T}, tapes::SparseVAFTape...) where {T}
    ops = AffineOperation.(tapes)
    A = hcat( (op.matrix for op in ops)...)
    b = +( (op.vector for op in ops)...)
    x = vcat( (tape.variables for tape in tapes)...)
    return SparseVAFTape([SparseAffineOperation(A, b)], x)
end

function MOIU.operate(::typeof(*), ::Type{T}, x::Number, tape::SparseVAFTape) where {T}
    d = MOI.output_dimension(tape)
    return add_operation!(tape, SparseAffineOperation(sparse(T(x)*I, d, d), Zero(d)))
end


# we do all pairs of `SparseVAFTape` and `AbstractVector{<:Number}`, and then do 3+ arguments by iterating
function MOIU.operate(::typeof(vcat), ::Type{T}, tape1::SparseVAFTape,
    tape2::SparseVAFTape) where {T}

    op1 = AffineOperation(tape1)
    op2 = AffineOperation(tape2)

    A = blockdiag(op1.matrix, op2.matrix)
    b = vcat(op1.vector, op2.vector)
    x = vcat(tape1.variables, tape2.variables)
    return SparseVAFTape([SparseAffineOperation(A, b)], x)
end

function MOIU.operate(::typeof(vcat), ::Type{T}, tape::SparseVAFTape,
    v::AbstractVector{<:Number}) where {T}
    op = AffineOperation(tape)
    n = length(v)
    m = size(op.matrix, 2) # bad for uniformscaling    
    Z = spzeros(T, n, m)
    A = vcat(op.matrix, Z)
    b = vcat(op.vector, v)
    return SparseVAFTape([SparseAffineOperation(A, b)], tape.variables)
end

function MOIU.operate(::typeof(vcat), ::Type{T}, v::AbstractVector{<:Number}, tape::SparseVAFTape) where {T}
    op = AffineOperation(tape)
    n = length(v)
    m = size(op.matrix, 2) # bad for uniformscaling    
    Z = spzeros(T, n, m)
    A = vcat(Z, op.matrix)
    b = vcat(v, op.vector)
    return SparseVAFTape([SparseAffineOperation(A, b)], tape.variables)
end


const SparseVAFTapeOrVec = Union{SparseVAFTape, AbstractVector{<:Number}}

function MOIU.operate(::typeof(vcat), ::Type{T}, arg1::SparseVAFTapeOrVec, arg2::SparseVAFTapeOrVec, arg3::SparseVAFTapeOrVec, args::Vararg{<:SparseVAFTapeOrVec}) where {T}
    all_args = (arg1, arg2, arg3, args...)
    foldl((a,b) -> MOIU.operate(vcat, T, a, b), all_args)
end
