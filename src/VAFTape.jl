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
            tape = collapse!(tape)
        end
        return tape
    end
end

function Base.isequal(a::VAFTape, b::VAFTape)
    isequal(a.variables, b.variables) && isequal(a.operations, b.operations)
end

function Base.hash(a::VAFTape{T}, h::UInt) where {T}
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



# # `MOIU.operate` methods for `VAFTape`


add_operation(tape::VAFTape, op::AffineOperation) = VAFTape((op, tape.operations...), tape.variables)

function MOIU.operate(::typeof(-), ::Type{T},
    tape::VAFTape) where {T}
    d = MOI.output_dimension(tape)
    return add_operation(tape, AffineOperation(-one(T)*I, Zero(d)))
end

# We need to revisit these... need to be able to handle any number of vectors, any number of `VAFTape`s, in any order....
###
function MOIU.operate(::typeof(+), ::Type{T}, v::AbstractVector,
                            tape::VAFTape) where {T}
    return add_operation(tape, AffineOperation(one(T)*I, v))
end

function MOIU.operate(::typeof(+), ::Type{T}, tape::VAFTape, vs::AbstractVector...) where {T}
    v = sum(vs)
    d = length(v)
    return add_operation(tape, AffineOperation(one(T)*I, v))
end

function MOIU.operate(::typeof(+), ::Type{T}, v::AbstractVector, tape::VAFTape, vs::AbstractVector...) where {T}
    return MOIU.operate(+, T, tape, v + sum(vs))
end

function MOIU.operate(::typeof(+), ::Type{T}, v::AbstractVector, tape::VAFTape) where {T}
    return MOIU.operate(+, T, tape, v)
end

function MOIU.operate(::typeof(+), ::Type{T}, tapes::VAFTape...) where {T}
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


function MOIU.operate(::typeof(-), ::Type{T}, tape::VAFTape,
    v::AbstractVector) where {T}
    return add_operation(tape, AffineOperation(one(T)*I, -v))
end

function MOIU.operate(::typeof(-), ::Type{T},
                    v::AbstractVector, tape::VAFTape) where {T}
    return add_operation(tape, AffineOperation(-one(T)*I, v))
end

function MOIU.operate(::typeof(*), ::Type{T}, A::AbstractMatrix,
    tape::VAFTape) where {T}
    return add_operation(tape, AffineOperation(A, Zero(size(A,1))))
end

function MOIU.operate(::typeof(sum), ::Type{T}, tape::VAFTape) where {T}
    d = MOI.output_dimension(tape)
    A = ones(T, 1, d) 
    return add_operation(tape, AffineOperation(A, Zero(size(A,1))))
end

function MOIU.operate(::typeof(vcat), ::Type{T}, tape1::VAFTape,
    tape2::VAFTape) where {T}
    op1 = AffineOperation(tape1)
    op2 = AffineOperation(tape2)
    M1 = op1.matrix
    M2 = op2.matrix
    Z1 = zeros(T, size(M2, 1), size(M1, 2))
    Z2 = zeros(T, size(M1, 1), size(M2, 2))
    A = [M1 Z1; Z2 M2]
    b = vcat(op1.vector, op2.vector)
    x = vcat(tape1.variables, tape2.variables)
    return VAFTape(tuple(AffineOperation(A, b)), x)
end

function MOIU.operate(::typeof(*), ::Type{T}, x::Number, tape::VAFTape) where {T}
    d = MOI.output_dimension(tape)
    return add_operation(tape, AffineOperation(T(x)*I, Zero(d)))
end
