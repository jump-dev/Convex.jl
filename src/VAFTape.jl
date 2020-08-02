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

# This should get upstreamed, or if it already is, version-gated
function LinearAlgebra.lmul!(a::Number, D::LinearAlgebra.Diagonal)
    lmul!(a, D.diag)
    return D
end

using MatrixChainMultiply

# How should `MatrixChainMultiply` treat the size of an `AffineOperation`?
# We will just take it to be the size of the matrix, disregarding sparsity etc.
MatrixChainMultiply.msize(a::AffineOperation) = size(a.matrix)

# If the matrix is a `UniformScaling` object that does not have a size, then resort to treating
# it as a square matrix the size of the vector.
MatrixChainMultiply.msize(a::AffineOperation{<:UniformScaling}) = (length(a.vector), length(a.vector))

# Override by `Convex.USE_CHAIN() = false`, which will trigger recompilation of the methods.
USE_CHAIN() = true

function AffineOperation!(tape::VAFTape)# -> AffineOperation
    if USE_CHAIN()
        op = matrixchainmultiply(compose!!, tape.operations...)
    else
        # left-to-right order
        t = tape.operations
        op, t = popfirst(t)
        while !isempty(t)
            newop, t = popfirst(t)
            op = compose!!(op, newop)
        end
    end
    return op
end

function compose!!(A::AffineOperation, B::AffineOperation)
    vec = try_mul!!(A.vector, A.matrix, B.vector, 1, 1)
    mat = try_mul!!(A.matrix, B.matrix)
    AffineOperation(mat, vec)
end

function collapse!(tape::VAFTape) # -> VAFTape
    op = AffineOperation!(tape)
    return VAFTape(tuple(op), tape.variables)
end



# What is the deal with this `try_mul!!`? The problem is that often
# one of the arguments is `LinearAlgebra.I`, or `Zero`, and we haven't actually
# allocated any memory there to do an inplace update. If we have allocated memory
# for the other argument, then we want to use that memory instead. The idea is that
# everything is scratch, and the return value just has to point to the right memory.
#
# We will need to be careful not to mutate user-inputted arrays!
"""
    try_mul!!(A, B)

Performs `A*B` and returns the result, possibly overwriting either `A` or `B`.
"""
function try_mul!!(A, B)
    rmul!(A, B)
end

LinearAlgebra.rmul!(z::Zero, a) = z
LinearAlgebra.lmul!(a, z::Zero) = z

"""
    try_mul!!(C, A, B, α, β)

Combined matrix-matrix or matrix-vector multiply-add A B α + C β, possibly `C` or `B`, returning the result.
"""
function try_mul!!(C, A, B, α, β)
    mul!(C, A, B, α, β)
end

function try_mul!!(C, A, ::Zero, α, β)
    try_mul!!(C, β)
end

function try_mul!!(A::UniformScaling, B)
    lmul!(A.λ, B)
end


function try_mul!!(A::Diagonal, B)
    lmul!(A, B)
end

function try_mul!!(A, B::UniformScaling)
    rmul!(A, B.λ)
end

function try_mul!!(A::UniformScaling, B::UniformScaling)
    (A.λ * B.λ)*I
end



# # `MOIU.operate` methods for `VAFTape`


add_operation(op::AffineOperation, tape::VAFTape) = VAFTape((op, tape.operations...), tape.variables)

function MOIU.operate(::typeof(-), ::Type{T},
    tape::VAFTape) where {T}
    d = MOI.output_dimension(tape)
    return add_operation(AffineOperation(-I, Zero(d)), tape)
end

function MOIU.operate(::typeof(+), ::Type{T}, v::AbstractVector,
                            tape::VAFTape) where {T}
    return add_operation(AffineOperation(I, v), tape)
end

function MOIU.operate(::typeof(-), ::Type{T}, tape::VAFTape,
    v::AbstractVector) where {T}
    return add_operation(AffineOperation(I, -v), tape)
end

function MOIU.operate(::typeof(-), ::Type{T},
                    v::AbstractVector, tape::VAFTape) where {T}
    return add_operation(AffineOperation(-I, v), tape)
end

function MOIU.operate(::typeof(*), ::Type{T}, A::AbstractMatrix,
    tape::VAFTape) where {T}
    return add_operation(AffineOperation(A, Zero(size(A,1))), tape)
end
