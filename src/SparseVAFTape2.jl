struct SparseAffineOperation2{T}
    matrix::SparseMatrixCSC{T}
    vector::Vector{T}
end


SparseAffineOperation2(A::AbstractSparseMatrix,b) = SparseAffineOperation2(A, b)
SparseAffineOperation2(A, b) = SparseAffineOperation2(sparse(A), b)


struct SparseVAFTape2{T}
    operations::Vector{SparseAffineOperation2{T}}
    variables::Vector{MOI.VariableIndex}
end

AffineOperation(op::SparseAffineOperation2) = AffineOperation(op.matrix, op.vector)
AffineOperation(op::SparseVAFTape2) = AffineOperation(SparseAffineOperation2(op))

const SparseVAFTape2OrVec = Union{SparseVAFTape2, AbstractVector{<:Real}}


MOI.output_dimension(v::SparseVAFTape2) = size(v.operations[1].matrix, 1)

function SparseAffineOperation2(tape::SparseVAFTape2)# -> AffineOperation
    return foldl(compose, tape.operations)
end

function compose(A::SparseAffineOperation2, B::SparseAffineOperation2)
    vec = A.vector + A.matrix * B.vector
    mat = A.matrix * B.matrix
    SparseAffineOperation2(mat, vec)
end

function collapse(sparse_tape::SparseVAFTape2) # -> SparseVAFTape2
    op = SparseAffineOperation2(sparse_tape)
    return SparseVAFTape2([op], sparse_tape.variables)
end
#### SparseVAFTape2

function add_operation(tape::SparseVAFTape2{T}, op::SparseAffineOperation2) where T
    tape2 = SparseVAFTape2(copy(tape.operations), tape.variables)
    pushfirst!(tape2.operations, op)
    return tape2::SparseVAFTape2{T}
end

function operate(::typeof(-), ::Type{T},
    tape::SparseVAFTape2) where {T <: Real}
    d = MOI.output_dimension(tape)
    return add_operation(tape, SparseAffineOperation2(sparse(-one(T)*I, d, d), zeros(T, d)))
end

function real_operate(::typeof(+), ::Type{T}, tape::SparseVAFTape2, v::AbstractVector{<:Real}) where {T <: Real}
    d = length(v)
    return add_operation(tape, SparseAffineOperation2(sparse(one(T)*I, d, d), v))
end

function real_operate(::typeof(+), ::Type{T}, v::AbstractVector{<:Real}, tape::SparseVAFTape2) where {T <: Real}
    return real_operate(+, T, tape, v)
end

function real_operate(::typeof(+), ::Type{T}, args::SparseVAFTape2OrVec...) where {T <: Real}
    vec_args = (a for a in args if a isa AbstractVector)
    tape_args = (a for a in args if a isa SparseVAFTape2)
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

function operate(::typeof(-), ::Type{T}, tape::SparseVAFTape2,
    v::AbstractVector{<:Real}) where {T <: Real}
    d = length(v)
    return add_operation(tape, SparseAffineOperation2(sparse(one(T)*I, d, d), -v))
end

function operate(::typeof(-), ::Type{T},
                    v::AbstractVector{<:Real}, tape::SparseVAFTape2) where {T <: Real}
    d = length(v)
    return add_operation(tape, SparseAffineOperation2(sparse(-one(T)*I, d, d), v))
end

function operate(::typeof(*), ::Type{T}, A::AbstractMatrix{<:Real},
    tape::SparseVAFTape2) where {T <: Real}
    return add_operation(tape, SparseAffineOperation2(A, zeros(T, size(A,1))))
end

function operate(::typeof(sum), ::Type{T}, tape::SparseVAFTape2) where {T <: Real}
    d = MOI.output_dimension(tape)
    # doesn't seem ideal for a sparse representation...
    A = ones(T, 1, d) 
    return add_operation(tape, SparseAffineOperation2(A, zeros(T, size(A,1))))
end


function operate(::typeof(+), ::Type{T}, tapes::SparseVAFTape2...) where {T <: Real}
    ops = AffineOperation.(tapes)
    A = hcat( (op.matrix for op in ops)...)
    b = +( (op.vector for op in ops)...)
    x = vcat( (tape.variables for tape in tapes)...)
    return SparseVAFTape2([SparseAffineOperation2(A, b)], x)
end

function operate(::typeof(*), ::Type{T}, x::Real, tape::SparseVAFTape2) where {T <: Real}
    d = MOI.output_dimension(tape)
    return add_operation(tape, SparseAffineOperation2(sparse(T(x)*I, d, d), zeros(T, d)))
end

# we do all pairs of `SparseVAFTape2` and `AbstractVector{<:Real}`, and then do 3+ arguments by iterating
function operate(::typeof(vcat), ::Type{T}, tape1::SparseVAFTape2,
    tape2::SparseVAFTape2) where {T <: Real}

    op1 = AffineOperation(tape1)
    op2 = AffineOperation(tape2)

    A = blockdiag(op1.matrix, op2.matrix)
    b = vcat(op1.vector, op2.vector)
    x = vcat(tape1.variables, tape2.variables)
    return SparseVAFTape2([SparseAffineOperation2(A, b)], x)
end

function operate(::typeof(vcat), ::Type{T}, tape::SparseVAFTape2,
    v::AbstractVector{<:Real}) where {T <: Real}
    op = AffineOperation(tape)
    n = length(v)
    m = size(op.matrix, 2) # bad for uniformscaling    
    Z = spzeros(T, n, m)
    A = vcat(op.matrix, Z)
    b = vcat(op.vector, v)
    return SparseVAFTape2([SparseAffineOperation2(A, b)], tape.variables)
end

function operate(::typeof(vcat), ::Type{T}, v::AbstractVector{<:Real}, tape::SparseVAFTape2) where {T <: Real}
    op = AffineOperation(tape)
    n = length(v)
    m = size(op.matrix, 2) # bad for uniformscaling    
    Z = spzeros(T, n, m)
    A = vcat(Z, op.matrix)
    b = vcat(v, op.vector)
    return SparseVAFTape2([SparseAffineOperation2(A, b)], tape.variables)
end


function operate(::typeof(vcat), ::Type{T}, arg1::SparseVAFTape2OrVec, arg2::SparseVAFTape2OrVec, arg3::SparseVAFTape2OrVec, args::Vararg{<:SparseVAFTape2OrVec}) where {T <: Real}
    all_args = (arg1, arg2, arg3, args...)
    foldl((a,b) -> operate(vcat, T, a, b), all_args)::SparseVAFTape2{T}
end


function operate(::typeof(+), ::Type{T}, tape1::SparseVAFTape2, tape2::SparseVAFTape2) where {T <: Real}
    @assert MOI.output_dimension(tape1) == MOI.output_dimension(tape2)
    op1 = AffineOperation(tape1)
    op2 = AffineOperation(tape2)

    if tape1.variables == tape2.variables
        op = SparseAffineOperation2(op1.matrix + op2.matrix, op1.vector + op2.vector)
        return SparseVAFTape2([op], tape1.variables)
    else
        mat = hcat(op1.matrix, op2.matrix)
        vec = op1.vector + op2.vector
        op = SparseAffineOperation2(mat,vec)
        return SparseVAFTape2([op], vcat(tape1.variables, tape2.variables))
    end
end
