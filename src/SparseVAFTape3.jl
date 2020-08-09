struct SparseAffineOperation3{T}
    matrix::Union{SparseMatrixCSC{T}, Diagonal{T, Vector{T}}}
    vector::Vector{T}
end


# SparseAffineOperation3(A::AbstractSparseMatrix,b) = SparseAffineOperation3(A, b)
# SparseAffineOperation3(A, b) = SparseAffineOperation3((A), b)


struct SparseVAFTape3{T}
    operations::Vector{SparseAffineOperation3{T}}
    variables::Vector{MOI.VariableIndex}
end

AffineOperation(op::SparseAffineOperation3) = AffineOperation(op.matrix, op.vector)
AffineOperation(op::SparseVAFTape3) = AffineOperation(SparseAffineOperation3(op))

const SparseVAFTape3OrVec = Union{SparseVAFTape3, AbstractVector{<:Real}}


MOI.output_dimension(v::SparseVAFTape3) = size(v.operations[1].matrix, 1)

function SparseAffineOperation3(tape::SparseVAFTape3)# -> AffineOperation
    return foldl(compose, tape.operations)
end

function compose(A::SparseAffineOperation3, B::SparseAffineOperation3)
    vec = A.vector + A.matrix * B.vector
    mat = A.matrix * B.matrix
    SparseAffineOperation3(mat, vec)
end

function collapse(sparse_tape::SparseVAFTape3) # -> SparseVAFTape3
    op = SparseAffineOperation3(sparse_tape)
    return SparseVAFTape3([op], sparse_tape.variables)
end
#### SparseVAFTape3

function add_operation(tape::SparseVAFTape3{T}, op::SparseAffineOperation3) where T
    tape2 = SparseVAFTape3(copy(tape.operations), tape.variables)
    pushfirst!(tape2.operations, op)
    return tape2::SparseVAFTape3{T}
end

function operate(::typeof(-), ::Type{T},
    tape::SparseVAFTape3) where {T <: Real}
    d = MOI.output_dimension(tape)
    return add_operation(tape, SparseAffineOperation3((-one(T)*I)(d), zeros(T, d)))
end

function real_operate(::typeof(+), ::Type{T}, tape::SparseVAFTape3, v::AbstractVector{<:Real}) where {T <: Real}
    d = length(v)
    return add_operation(tape, SparseAffineOperation3((one(T)*I)(d), v))
end

function real_operate(::typeof(+), ::Type{T}, v::AbstractVector{<:Real}, tape::SparseVAFTape3) where {T <: Real}
    return real_operate(+, T, tape, v)
end

function real_operate(::typeof(+), ::Type{T}, args::SparseVAFTape3OrVec...) where {T <: Real}
    vec_args = (a for a in args if a isa AbstractVector)
    tape_args = (a for a in args if a isa SparseVAFTape3)
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

function operate(::typeof(-), ::Type{T}, tape::SparseVAFTape3,
    v::AbstractVector{<:Real}) where {T <: Real}
    d = length(v)
    return add_operation(tape, SparseAffineOperation3((one(T)*I)(d), -v))
end

function operate(::typeof(-), ::Type{T},
                    v::AbstractVector{<:Real}, tape::SparseVAFTape3) where {T <: Real}
    d = length(v)
    return add_operation(tape, SparseAffineOperation3((-one(T)*I)(d), v))
end

function operate(::typeof(*), ::Type{T}, A::AbstractMatrix{<:Real},
    tape::SparseVAFTape3) where {T <: Real}
    return add_operation(tape, SparseAffineOperation3(A, zeros(T, size(A,1))))
end

function operate(::typeof(sum), ::Type{T}, tape::SparseVAFTape3) where {T <: Real}
    d = MOI.output_dimension(tape)
    # doesn't seem ideal for a sparse representation?
    A = sparse(ones(T, 1, d))
    return add_operation(tape, SparseAffineOperation3(A, zeros(T, size(A,1))))
end


function operate(::typeof(+), ::Type{T}, tapes::SparseVAFTape3...) where {T <: Real}
    ops = AffineOperation.(tapes)
    A = hcat( (op.matrix for op in ops)...)
    b = +( (op.vector for op in ops)...)
    x = vcat( (tape.variables for tape in tapes)...)
    return SparseVAFTape3([SparseAffineOperation3(A, b)], x)
end

function operate(::typeof(*), ::Type{T}, x::Real, tape::SparseVAFTape3) where {T <: Real}
    d = MOI.output_dimension(tape)
    return add_operation(tape, SparseAffineOperation3((T(x)*I)(d), zeros(T, d)))
end

_sparse(mat) = mat isa Diagonal ? sparse(mat) : mat

# we do all pairs of `SparseVAFTape3` and `AbstractVector{<:Real}`, and then do 3+ arguments by iterating
function operate(::typeof(vcat), ::Type{T}, tape1::SparseVAFTape3,
    tape2::SparseVAFTape3) where {T <: Real}

    op1 = AffineOperation(tape1)
    op2 = AffineOperation(tape2)

    A = blockdiag(_sparse(op1.matrix), _sparse(op2.matrix))
    b = vcat(op1.vector, op2.vector)
    x = vcat(tape1.variables, tape2.variables)
    return SparseVAFTape3([SparseAffineOperation3(A, b)], x)
end

function operate(::typeof(vcat), ::Type{T}, tape::SparseVAFTape3,
    v::AbstractVector{<:Real}) where {T <: Real}
    op = AffineOperation(tape)
    n = length(v)
    m = size(op.matrix, 2) # bad for uniformscaling    
    Z = spzeros(T, n, m)
    A = vcat(op.matrix, Z)
    b = vcat(op.vector, v)
    return SparseVAFTape3([SparseAffineOperation3(A, b)], tape.variables)
end

function operate(::typeof(vcat), ::Type{T}, v::AbstractVector{<:Real}, tape::SparseVAFTape3) where {T <: Real}
    op = AffineOperation(tape)
    n = length(v)
    m = size(op.matrix, 2) # bad for uniformscaling    
    Z = spzeros(T, n, m)
    A = vcat(Z, op.matrix)
    b = vcat(v, op.vector)
    return SparseVAFTape3([SparseAffineOperation3(A, b)], tape.variables)
end


function operate(::typeof(vcat), ::Type{T}, arg1::SparseVAFTape3OrVec, arg2::SparseVAFTape3OrVec, arg3::SparseVAFTape3OrVec, args::Vararg{<:SparseVAFTape3OrVec}) where {T <: Real}
    all_args = (arg1, arg2, arg3, args...)
    foldl((a,b) -> operate(vcat, T, a, b), all_args)::SparseVAFTape3{T}
end


function operate(::typeof(+), ::Type{T}, tape1::SparseVAFTape3, tape2::SparseVAFTape3) where {T <: Real}
    @assert MOI.output_dimension(tape1) == MOI.output_dimension(tape2)
    op1 = AffineOperation(tape1)
    op2 = AffineOperation(tape2)

    if tape1.variables == tape2.variables
        op = SparseAffineOperation3(op1.matrix + op2.matrix, op1.vector + op2.vector)
        return SparseVAFTape3([op], tape1.variables)
    else
        mat = hcat(op1.matrix, op2.matrix)
        vec = op1.vector + op2.vector
        op = SparseAffineOperation3(mat,vec)
        return SparseVAFTape3([op], vcat(tape1.variables, tape2.variables))
    end
end
