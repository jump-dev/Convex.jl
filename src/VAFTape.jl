# My original implementation of the VAF tape idea.
# `AffineOperation` is parametrized on both the `matrix` type and `vector` type
# allowing use of `LinearAlgebra.I` and special vector types, etc.
# This is nice because we can use optimized implementations when relevant,
# but can be more cumbersome / error prone to ensure optimized dispatches are all hit,
# and because type stability can be an issue.

struct AffineOperation{M, V}
    matrix::M
    vector::V
end


function Base.isequal(A::AffineOperation, B::AffineOperation)
    isequal(A.vector, B.vector) && isequal(A.matrix, B.matrix)
end

function Base.hash(a::AffineOperation{M, V}, h::UInt) where {M, V}
    hash(a.matrix, hash(a.vector, hash(:AffineOperation, h)))
end

# configurable parameter which is a compile time constant
COLLAPSE_DEPTH() = 15

struct VAFTape{T <: Tuple}
    operations::T
    variables::Vector{MOI.VariableIndex}

    function VAFTape(operations::T, variables::Vector{MOI.VariableIndex}) where {T <: Tuple}
        tape = new{T}(operations, variables)
        # don't let our tuples get too long!
        if length(operations) > COLLAPSE_DEPTH()
            tape = collapse(tape)
        end
        return tape
    end
end
const VAFTapeOrVec = Union{VAFTape, AbstractVector{<:Real}}

function Base.isequal(a::VAFTape, b::VAFTape)
    isequal(a.variables, b.variables) && isequal(a.operations, b.operations)
end

function Base.hash(a::VAFTape{T}, h::UInt) where {T <: Real}
    hash(a.operations, hash(a.variables, hash(:VAFTape, h)))
end


function MOI.output_dimension(v::VAFTape)
    t = v.operations
    while !isempty(t)
        op, t = popfirst(t)
        if !(op.matrix isa UniformScaling)
            return size(op.matrix, 1)
        end
    end
    return length(v.variables)
end


function popfirst(t::Tuple)
    return first(t), Base.tail(t)
end


using MatrixChainMultiply

# How should `MatrixChainMultiply` treat the size of an `AffineOperation`?
# We will just take it to be the size of the matrix, disregarding sparsity etc.
MatrixChainMultiply.msize(a::AffineOperation) = size(a.matrix)

# If the matrix is a `UniformScaling` object that does not have a size, then resort to treating
# it as a square matrix the size of the vector.
MatrixChainMultiply.msize(a::AffineOperation{<:UniformScaling}) = (length(a.vector), length(a.vector))

# Override by `Convex.USE_CHAIN() = false`, which will trigger recompilation of the methods.
USE_CHAIN() = false

function AffineOperation(tape::VAFTape)# -> AffineOperation
    if USE_CHAIN()
        # this will be inherently type-unstable...
        op = matrixchainmultiply(compose, tape.operations...)
    else
        # left-to-right order; should be type stable
        t = tape.operations
        op, t = popfirst(t)
        while !isempty(t)
            newop, t = popfirst(t)
            op = compose(op, newop)
        end
    end
    return op
end


function compose(A::AffineOperation, B::AffineOperation)
    vec = A.vector + A.matrix * B.vector
    mat = A.matrix * B.matrix
    AffineOperation(mat, vec)
end

function collapse(tape::VAFTape) # -> VAFTape
    op = AffineOperation(tape)
    return VAFTape(tuple(op), tape.variables)
end



# # `operate` methods for `VAFTape`


add_operation(tape::VAFTape, op::AffineOperation) = VAFTape((op, tape.operations...), tape.variables)

function operate(::typeof(-), ::Type{T},
    tape::VAFTape) where {T <: Real}
    d = MOI.output_dimension(tape)
    return add_operation(tape, AffineOperation(-one(T)*I, zeros(T, d)))
end

# We need to revisit these... need to be able to handle any number of vectors, any number of `VAFTape`s, in any order....
###

function real_operate(::typeof(+), ::Type{T}, v::AbstractVector,
                            tape::VAFTape) where {T <: Real}
    return add_operation(tape, AffineOperation(one(T)*I, v))
end

function real_operate(::typeof(+), ::Type{T}, tape::VAFTape, v::AbstractVector{T}) where {T <: Real}
    d = length(v)
    return add_operation(tape, AffineOperation(one(T)*I, v))
end

# function real_operate(::typeof(+), ::Type{T}, v::AbstractVector, tape::VAFTape, vs::AbstractVector...) where {T <: Real}
    # return real_operate(+, T, tape, v + sum(vs))
# end


function real_operate(::typeof(+), ::Type{T}, args::VAFTapeOrVec...) where {T <: Real}
    vec_args = (a for a in args if a isa AbstractVector)
    tape_args = (a for a in args if a isa VAFTape)
    if isempty(tape_args)
        return sum(vec_args)
    else
        tape = foldl((a,b) -> real_operate(+, T, a, b), tape_args)
    end
    if isempty(vec_args)
        return tape
    else
        v = sum(vec_args)
        return real_operate(+, T, tape, v)
    end
end


function real_operate(::typeof(+), ::Type{T}, tapes::VAFTape...) where {T <: Real}
    ops = AffineOperation.(tapes)
    if all(op -> op.matrix isa UniformScaling, ops)
        # if they are all UniformScaling, we can't `hcat` since it doesn't know the size
        # but we know the size, so we'll instantiate an `I(d)` and use that
        d= length(ops[end].vector)
        λ = ops[end].matrix.λ 
        A = hcat( (op.matrix for op in ops[1:end-1])..., λ*I(d))
    else
        A = hcat( (op.matrix for op in ops)...)
    end
    b = +( (op.vector for op in ops)...)
    x = vcat( (tape.variables for tape in tapes)...)
    return VAFTape(tuple(AffineOperation(A, b)), x)
end

####


function operate(::typeof(-), ::Type{T}, tape::VAFTape,
    v::AbstractVector) where {T <: Real}
    return add_operation(tape, AffineOperation(one(T)*I, -v))
end

function operate(::typeof(-), ::Type{T},
                    v::AbstractVector, tape::VAFTape) where {T <: Real}
    return add_operation(tape, AffineOperation(-one(T)*I, v))
end

function operate(::typeof(*), ::Type{T}, A::AbstractMatrix,
    tape::VAFTape) where {T <: Real}
    return add_operation(tape, AffineOperation(A, zeros(T, size(A,1))))
end

function operate(::typeof(sum), ::Type{T}, tape::VAFTape) where {T <: Real}
    d = MOI.output_dimension(tape)
    A = ones(T, 1, d) 
    return add_operation(tape, AffineOperation(A, zeros(T, size(A,1))))
end


# we do all pairs of `SparseVAFTape` and `AbstractVector{<:Real}`, and then do 3+ arguments by iterating


function operate(::typeof(vcat), ::Type{T}, tape::VAFTape,
    v::AbstractVector{<:Real}) where {T <: Real}
    op = AffineOperation(tape)
    n = length(v)
    if op.matrix isa UniformScaling
        m = n
    else
        m = size(op.matrix, 2)  
    end
    Z = spzeros(T, n, m)
    A = vcat(op.matrix, Z)
    b = vcat(op.vector, v)
    return VAFTape(tuple(AffineOperation(A, b)), tape.variables)
end

function operate(::typeof(vcat), ::Type{T}, v::AbstractVector{<:Real}, tape::VAFTape) where {T <: Real}
    op = AffineOperation(tape)
    n = length(v)
    if op.matrix isa UniformScaling
        m = n
    else
        m = size(op.matrix, 2)  
    end
    Z = spzeros(T, n, m)
    A = vcat(Z, op.matrix)
    b = vcat(v, op.vector)
    return VAFTape(tuple(AffineOperation(A, b)), tape.variables)
end


function operate(::typeof(vcat), ::Type{T}, arg1::VAFTapeOrVec, arg2::VAFTapeOrVec, arg3::VAFTapeOrVec, args::Vararg{<:VAFTapeOrVec}) where {T <: Real}
    all_args = (arg1, arg2, arg3, args...)
    foldl((a,b) -> operate(vcat, T, a, b), all_args)
end

function operate(::typeof(vcat), ::Type{T}, tape1::VAFTape, tape2::VAFTape) where {T <: Real}
    op1 = AffineOperation(tape1)
    op2 = AffineOperation(tape2)
    M1 = op1.matrix
    M2 = op2.matrix
    if M2 isa UniformScaling
        m2 = n2 = length(op2.vector)
    else
        m2, n2 = size(M2)
    end
    if M1 isa UniformScaling
        m1 = n1 = length(op1.vector)
    else
        m1, n1 = size(M1)
    end
    Z1 = zeros(T, m1, n2)
    Z2 = zeros(T, m2, n1)
    A = [M1 Z1; Z2 M2]
    b = vcat(op1.vector, op2.vector)
    x = vcat(tape1.variables, tape2.variables)
    return VAFTape(tuple(AffineOperation(A, b)), x)
end

function operate(::typeof(*), ::Type{T}, x::Real, tape::VAFTape) where {T <: Real}
    d = MOI.output_dimension(tape)
    return add_operation(tape, AffineOperation(T(x)*I, zeros(T, d)))
end
