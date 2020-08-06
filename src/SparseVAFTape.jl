# An alternate implementation of the VAF tape concept:
# Here, we represent each affine transformation by a single sparse matrix
# of uniform type. This helps avoid type issues, since every operation has the same type,
# and does not rely on correctly hitting optimized dispatches.
# On the other hand, we can't use optimized dispatches to help speed things up when
# we have extra structural information.

struct SparseAffineOperation{T <: AbstractSparseMatrix{<:Real}}
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

const SparseVAFTapeOrVec = Union{SparseVAFTape, AbstractVector{<:Real}}


MOI.output_dimension(v::SparseVAFTape) = size(v.operations[1].matrix, 1) - 1


_unwrap(a::SparseAffineOperation) = a.matrix
_unwrap(a::AbstractSparseArray) = a

function AffineOperation(sparse_tape::SparseVAFTape)
    sparse_tape = collapse(sparse_tape)
    mat = only(sparse_tape.operations).matrix
    b = mat[1:end-1, end]
    A = mat[1:end-1, 1:end-1]
    AffineOperation(A, b)
end

function SparseAffineOperation(sparse_tape::SparseVAFTape)
    sparse_tape = collapse(sparse_tape)
    return sparse_tape.operations[]
end

function collapse(sparse_tape::SparseVAFTape) # -> SparseVAFTape
    # @show size.(_unwrap.(sparse_tape.operations))
    if length(sparse_tape.operations) > 1
        mat = foldl((a,b) -> _unwrap(a) * _unwrap(b), sparse_tape.operations)
    else
        mat = only(sparse_tape.operations).matrix
    end
    op = SparseAffineOperation(mat)
    return SparseVAFTape([op], sparse_tape.variables)
end
#### SparseVAFTape

function add_operation(tape::SparseVAFTape{T}, op::SparseAffineOperation) where T
    tape2 = SparseVAFTape(copy(tape.operations), tape.variables)
    pushfirst!(tape2.operations, op)
    return tape2::SparseVAFTape{T}
end

function operate(::typeof(-), ::Type{T},
    tape::SparseVAFTape) where {T <: Real}
    d = MOI.output_dimension(tape)
    return add_operation(tape, SparseAffineOperation(sparse(-one(T)*I, d, d), zeros(T, d)))
end

function real_operate(::typeof(+), ::Type{T}, tape::SparseVAFTape, v::AbstractVector{<:Real}) where {T <: Real}
    d = length(v)
    return add_operation(tape, SparseAffineOperation(sparse(one(T)*I, d, d), v))
end

function real_operate(::typeof(+), ::Type{T}, v::AbstractVector{<:Real}, tape::SparseVAFTape) where {T <: Real}
    return real_operate(+, T, tape, v)
end

function real_operate(::typeof(+), ::Type{T}, args::SparseVAFTapeOrVec...) where {T <: Real}
    vec_args = (a for a in args if a isa AbstractVector)
    tape_args = (a for a in args if a isa SparseVAFTape)
    if isempty(tape_args)
        error()
    else
        tape = foldl((a,b) -> operate(+, T, a, b), tape_args)
    end
    if isempty(vec_args)
        return tape
    else
        v = sum(vec_args)
        return operate(+, T, tape, v)
    end
end

function operate(::typeof(-), ::Type{T}, tape::SparseVAFTape,
    v::AbstractVector{<:Real}) where {T <: Real}
    d = length(v)
    return add_operation(tape, SparseAffineOperation(sparse(one(T)*I, d, d), -v))
end

function operate(::typeof(-), ::Type{T},
                    v::AbstractVector{<:Real}, tape::SparseVAFTape) where {T <: Real}
    d = length(v)
    return add_operation(tape, SparseAffineOperation(sparse(-one(T)*I, d, d), v))
end

function operate(::typeof(*), ::Type{T}, A::AbstractMatrix{<:Real},
    tape::SparseVAFTape) where {T <: Real}
    return add_operation(tape, SparseAffineOperation(A, zeros(T, size(A,1))))
end

function operate(::typeof(sum), ::Type{T}, tape::SparseVAFTape) where {T <: Real}
    d = MOI.output_dimension(tape)
    # doesn't seem ideal for a sparse representation...
    A = ones(T, 1, d) 
    return add_operation(tape, SparseAffineOperation(A, zeros(T, size(A,1))))
end


function operate(::typeof(+), ::Type{T}, tapes::SparseVAFTape...) where {T <: Real}
    ops = AffineOperation.(tapes)
    A = hcat( (op.matrix for op in ops)...)
    b = +( (op.vector for op in ops)...)
    x = vcat( (tape.variables for tape in tapes)...)
    return SparseVAFTape([SparseAffineOperation(A, b)], x)
end

function operate(::typeof(*), ::Type{T}, x::Real, tape::SparseVAFTape) where {T <: Real}
    d = MOI.output_dimension(tape)
    return add_operation(tape, SparseAffineOperation(sparse(T(x)*I, d, d), zeros(T, d)))
end

# we do all pairs of `SparseVAFTape` and `AbstractVector{<:Real}`, and then do 3+ arguments by iterating
function operate(::typeof(vcat), ::Type{T}, tape1::SparseVAFTape,
    tape2::SparseVAFTape) where {T <: Real}

    op1 = AffineOperation(tape1)
    op2 = AffineOperation(tape2)

    A = blockdiag(op1.matrix, op2.matrix)
    b = vcat(op1.vector, op2.vector)
    x = vcat(tape1.variables, tape2.variables)
    return SparseVAFTape([SparseAffineOperation(A, b)], x)
end

function operate(::typeof(vcat), ::Type{T}, tape::SparseVAFTape,
    v::AbstractVector{<:Real}) where {T <: Real}
    op = AffineOperation(tape)
    n = length(v)
    m = size(op.matrix, 2) # bad for uniformscaling    
    Z = spzeros(T, n, m)
    A = vcat(op.matrix, Z)
    b = vcat(op.vector, v)
    return SparseVAFTape([SparseAffineOperation(A, b)], tape.variables)
end

function operate(::typeof(vcat), ::Type{T}, v::AbstractVector{<:Real}, tape::SparseVAFTape) where {T <: Real}
    op = AffineOperation(tape)
    n = length(v)
    m = size(op.matrix, 2) # bad for uniformscaling    
    Z = spzeros(T, n, m)
    A = vcat(Z, op.matrix)
    b = vcat(v, op.vector)
    return SparseVAFTape([SparseAffineOperation(A, b)], tape.variables)
end


function operate(::typeof(vcat), ::Type{T}, arg1::SparseVAFTapeOrVec, arg2::SparseVAFTapeOrVec, arg3::SparseVAFTapeOrVec, args::Vararg{<:SparseVAFTapeOrVec}) where {T <: Real}
    all_args = (arg1, arg2, arg3, args...)
    foldl((a,b) -> operate(vcat, T, a, b), all_args)::SparseVAFTape{T}
end


function operate(::typeof(+), ::Type{T}, tape1::SparseVAFTape, tape2::SparseVAFTape) where {T <: Real}
    @assert MOI.output_dimension(tape1) == MOI.output_dimension(tape2)
    op1 = AffineOperation(tape1)
    op2 = AffineOperation(tape2)

    if tape1.variables == tape2.variables
        op = SparseAffineOperation(op1.matrix + op2.matrix, op1.vector + op2.vector)
        return SparseVAFTape([op], tape1.variables)
    else
        mat = hcat(op1.matrix, op2.matrix)
        vec = op1.vector + op2.vector
        op = SparseAffineOperation(mat,vec)
        return SparseVAFTape([op], vcat(tape1.variables, tape2.variables))
    end
end
