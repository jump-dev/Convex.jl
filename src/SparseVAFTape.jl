struct SparseAffineOperation{T}
    matrix::SparseMatrixCSC{T}
    vector::Vector{T}
end

function SparseAffineOperation(A::AbstractSparseMatrix, b)
    return SparseAffineOperation(SparseMatrixCSC(A), b)
end
SparseAffineOperation(A, b) = SparseAffineOperation(sparse(A), b)

struct SparseVAFTape{T}
    operations::Vector{SparseAffineOperation{T}}
    variables::Vector{MOI.VariableIndex}
end

function AffineOperation(op::SparseAffineOperation)
    return AffineOperation(op.matrix, op.vector)
end
function AffineOperation(op::SparseVAFTape)
    return AffineOperation(SparseAffineOperation(op))
end

const SparseVAFTapeOrVec = Union{SparseVAFTape,AbstractVector{<:Real}}

MOI.output_dimension(v::SparseVAFTape) = size(v.operations[1].matrix, 1)

function SparseAffineOperation(tape::SparseVAFTape)# -> AffineOperation
    return foldl(compose, tape.operations)
end

function compose(A::SparseAffineOperation, B::SparseAffineOperation)
    vec = A.vector + A.matrix * B.vector
    mat = A.matrix * B.matrix
    return SparseAffineOperation(mat, vec)
end

function collapse(sparse_tape::SparseVAFTape) # -> SparseVAFTape
    op = SparseAffineOperation(sparse_tape)
    return SparseVAFTape([op], sparse_tape.variables)
end
#### SparseVAFTape

function add_operation(
    tape::SparseVAFTape{T},
    op::SparseAffineOperation,
) where {T}
    tape2 = SparseVAFTape(copy(tape.operations), tape.variables)
    pushfirst!(tape2.operations, op)
    return tape2::SparseVAFTape{T}
end

function operate(::typeof(-), ::Type{T}, tape::SparseVAFTape) where {T<:Real}
    d = MOI.output_dimension(tape)
    return add_operation(
        tape,
        SparseAffineOperation(sparse(-one(T) * I, d, d), zeros(T, d)),
    )
end

function real_operate(
    ::typeof(+),
    ::Type{T},
    tape::SparseVAFTape,
    v::AbstractVector{<:Real},
) where {T<:Real}
    d = length(v)
    return add_operation(
        tape,
        SparseAffineOperation(sparse(one(T) * I, d, d), v),
    )
end

function real_operate(
    ::typeof(+),
    ::Type{T},
    v::AbstractVector{<:Real},
    tape::SparseVAFTape,
) where {T<:Real}
    return real_operate(+, T, tape, v)
end

function real_operate(
    ::typeof(+),
    ::Type{T},
    args::SparseVAFTapeOrVec...,
) where {T<:Real}
    vec_args = (a for a in args if a isa AbstractVector)
    tape_args = (a for a in args if a isa SparseVAFTape)
    if isempty(tape_args)
        error()
    else
        tape = foldl((a, b) -> operate(+, T, a, b), tape_args)
    end
    if isempty(vec_args)
        return tape
    else
        v = sum(vec_args)
        return operate(+, T, tape, v)
    end
end

function operate(
    ::typeof(-),
    ::Type{T},
    tape::SparseVAFTape,
    v::AbstractVector{<:Real},
) where {T<:Real}
    d = length(v)
    return add_operation(
        tape,
        SparseAffineOperation(sparse(one(T) * I, d, d), -v),
    )
end

function operate(
    ::typeof(-),
    ::Type{T},
    v::AbstractVector{<:Real},
    tape::SparseVAFTape,
) where {T<:Real}
    d = length(v)
    return add_operation(
        tape,
        SparseAffineOperation(sparse(-one(T) * I, d, d), v),
    )
end

function operate(
    ::typeof(*),
    ::Type{T},
    A::AbstractMatrix{<:Real},
    tape::SparseVAFTape,
) where {T<:Real}
    return add_operation(tape, SparseAffineOperation(A, zeros(T, size(A, 1))))
end

function operate(::typeof(sum), ::Type{T}, tape::SparseVAFTape) where {T<:Real}
    d = MOI.output_dimension(tape)
    # doesn't seem ideal for a sparse representation...
    A = ones(T, 1, d)
    return add_operation(tape, SparseAffineOperation(A, zeros(T, size(A, 1))))
end

function operate(
    ::typeof(+),
    ::Type{T},
    tapes::SparseVAFTape...,
) where {T<:Real}
    ops = AffineOperation.(tapes)
    A = hcat((op.matrix for op in ops)...)
    b = +((op.vector for op in ops)...)
    x = vcat((tape.variables for tape in tapes)...)
    return SparseVAFTape([SparseAffineOperation(A, b)], x)
end

function operate(
    ::typeof(*),
    ::Type{T},
    x::Real,
    tape::SparseVAFTape,
) where {T<:Real}
    d = MOI.output_dimension(tape)
    return add_operation(
        tape,
        SparseAffineOperation(sparse(T(x) * I, d, d), zeros(T, d)),
    )
end

# we do all pairs of `SparseVAFTape` and `AbstractVector{<:Real}`, and then do 3+ arguments by iterating
function operate(
    ::typeof(vcat),
    ::Type{T},
    tape1::SparseVAFTape,
    tape2::SparseVAFTape,
) where {T<:Real}
    op1 = AffineOperation(tape1)
    op2 = AffineOperation(tape2)

    A = blockdiag(op1.matrix, op2.matrix)
    b = vcat(op1.vector, op2.vector)
    x = vcat(tape1.variables, tape2.variables)
    return SparseVAFTape([SparseAffineOperation(A, b)], x)
end

function operate(
    ::typeof(vcat),
    ::Type{T},
    tape::SparseVAFTape,
    v::AbstractVector{<:Real},
) where {T<:Real}
    op = AffineOperation(tape)
    n = length(v)
    m = size(op.matrix, 2) # bad for uniformscaling
    Z = spzeros(T, n, m)
    A = vcat(op.matrix, Z)
    b = vcat(op.vector, v)
    return SparseVAFTape([SparseAffineOperation(A, b)], tape.variables)
end

function operate(
    ::typeof(vcat),
    ::Type{T},
    v::AbstractVector{<:Real},
    tape::SparseVAFTape,
) where {T<:Real}
    op = AffineOperation(tape)
    n = length(v)
    m = size(op.matrix, 2) # bad for uniformscaling
    Z = spzeros(T, n, m)
    A = vcat(Z, op.matrix)
    b = vcat(v, op.vector)
    return SparseVAFTape([SparseAffineOperation(A, b)], tape.variables)
end

function operate(
    ::typeof(vcat),
    ::Type{T},
    arg1::SparseVAFTapeOrVec,
    arg2::SparseVAFTapeOrVec,
    arg3::SparseVAFTapeOrVec,
    args::Vararg{<:SparseVAFTapeOrVec},
) where {T<:Real}
    all_args = (arg1, arg2, arg3, args...)
    return foldl((a, b) -> operate(vcat, T, a, b), all_args)::SparseVAFTape{T}
end

function operate(
    ::typeof(+),
    ::Type{T},
    tape1::SparseVAFTape,
    tape2::SparseVAFTape,
) where {T<:Real}
    @assert MOI.output_dimension(tape1) == MOI.output_dimension(tape2)
    op1 = AffineOperation(tape1)
    op2 = AffineOperation(tape2)

    if tape1.variables == tape2.variables
        op = SparseAffineOperation(
            op1.matrix + op2.matrix,
            op1.vector + op2.vector,
        )
        return SparseVAFTape([op], tape1.variables)
    else
        mat = hcat(op1.matrix, op2.matrix)
        vec = op1.vector + op2.vector
        op = SparseAffineOperation(mat, vec)
        return SparseVAFTape([op], vcat(tape1.variables, tape2.variables))
    end
end
