# configurable parameter which is a compile time constant
COLLAPSE_DEPTH() = 15

struct AffineOperation{M, V}
    matrix::M
    vector::V
end

struct SparseAffineOperation{T <: AbstractSparseMatrix}
    matrix::T
end

using SuiteSparse
USE_CHOLMOD() = false


SparseAffineOperation(A,b::AbstractSparseVector) = SparseAffineOperation(sparse(A), b)
SparseAffineOperation(A::AbstractSparseMatrix,b) = SparseAffineOperation(A, sparse(b))
SparseAffineOperation(A, b) = SparseAffineOperation(sparse(A), sparse(b))

function SparseAffineOperation(A::Union{SparseMatrixCSC, SuiteSparse.CHOLMOD.Sparse}, b::Union{SparseVector, SuiteSparse.CHOLMOD.Sparse})
    T = eltype(A)
    # construct a sparse matrix representation of `Ax+b`
    # via `[A b; 0 1] * [x; 1] == [Ax+b; 1]`.
    n, m = size(A)
    A = sparse(A)
    mat = [A b; transpose(spzeros(T, m)) one(T)]
    if USE_CHOLMOD()
        mat = SuiteSparse.CHOLMOD.Sparse(mat)
    end
    SparseAffineOperation(mat)
end

struct SparseVAFTape{T}
    operations::Vector{SparseAffineOperation{T}}
    variables::Vector{MOI.VariableIndex}
end

_unwrap(a::SparseAffineOperation) = a.matrix
_unwrap(a::AbstractSparseArray) = a

function AffineOperation!(sparse_tape::SparseVAFTape)
    mat = foldl((a,b) -> _unwrap(a) * _unwrap(b), sparse_tape.operations)
    b = mat[1:end-1, end]
    A = mat[1:end-1, 1:end-1]
    AffineOperation(A, b)
end


function Base.isequal(A::AffineOperation, B::AffineOperation)
    isequal(A.vector, B.vector) && isequal(A.matrix, B.matrix)
end

function Base.hash(a::AffineOperation{M, V}, h::UInt) where {M, V}
    hash(a.matrix, hash(a.vector, hash(:AffineOperation, h)))
end


# This is a variant of MOI.VectorAffineFunction which represents
# the transformation `matrix * variables + vector` lazily.
struct VectorAffineFunctionAsMatrix{M,B,V}
    aff::AffineOperation{M, B}
    variables::V
end

function Base.isequal(A::VectorAffineFunctionAsMatrix, B::VectorAffineFunctionAsMatrix)
    isequal(A.variables, B.variables) && isequal(A.aff, B.aff)
end


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

const VAFTapes = Union{VAFTape, SparseVAFTape}


function Base.isequal(a::VAFTape, b::VAFTape)
    isequal(a.variables, b.variables) && isequal(a.operations, b.operations)
end

function Base.hash(a::VAFTape{T}, h::UInt) where {T}
    hash(a.operations, hash(a.variables, hash(:VAFTape, h)))
end


# A simple type representing a vector of zeros. Maybe should include the size or use FillArrays or similar.
struct Zero
    len::Int
 end

Base.:(+)(a, ::Zero) = a
Base.:(+)(::Zero, a) = a
Base.:(-)(::Zero, a) = -a
Base.:(-)(a, ::Zero) = a
Base.:(-)(z::Zero) = z
Base.:(*)(A, z::Zero) = Zero(size(A, 1))
Base.size(z::Zero) = (z.len,)
Base.length(z::Zero) = z.len
SparseArrays.sparse(z::Zero) = spzeros(z.len)

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

MOI.output_dimension(v::SparseVAFTape) = size(v.operations[1].matrix, 1) - 1

function MOI.output_dimension(v::VectorAffineFunctionAsMatrix)
    return size(v.matrix, 1)
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

function to_vaf(tape::VAFTapes)
    op = AffineOperation!(tape)
    to_vaf(VectorAffineFunctionAsMatrix(op, tape.variables))
end


# convert to a usual VAF
function to_vaf(vaf_as_matrix::VectorAffineFunctionAsMatrix{<:AbstractSparseArray})
    T = eltype(vaf_as_matrix.aff.matrix)
    I, J, V = findnz(vaf_as_matrix.aff.matrix)
    vats = MOI.VectorAffineTerm{T}[]
    for n in eachindex(I, J, V)
        i = I[n]
        j = J[n]
        v = V[n]
        push!(vats,
              MOI.VectorAffineTerm{T}(i,
                                      MOI.ScalarAffineTerm{T}(v,
                                                              vaf_as_matrix.variables[j])))
    end

    return MOI.VectorAffineFunction{T}(vats, vaf_as_matrix.aff.vector)
end

function to_vaf(vaf_as_matrix::VectorAffineFunctionAsMatrix{<:AbstractMatrix})
    T = eltype(vaf_as_matrix.aff.matrix)
    vats = MOI.VectorAffineTerm{T}[]
    M = vaf_as_matrix.aff.matrix
    for i in 1:size(M, 1)
        for j in 1:size(M, 2)
            iszero(M[i, j]) && continue
            push!(vats,
                  MOI.VectorAffineTerm{T}(i,
                                          MOI.ScalarAffineTerm{T}(M[i, j],
                                                                  vaf_as_matrix.variables[j])))
        end
    end
    return MOI.VectorAffineFunction{T}(vats, vaf_as_matrix.aff.vector)
end

# method for adding constraints and coverting to standard VAFs as needed
MOI_add_constraint(model, f, set) = MOI.add_constraint(model, f, set)

function MOI_add_constraint(model, f::VectorAffineFunctionAsMatrix, set)
    return MOI.add_constraint(model, to_vaf(f), set)
end

function MOI_add_constraint(model, f::VAFTapes, set)
    return MOI.add_constraint(model, to_vaf(f), set)
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

#### SparseVAFTape

function add_operation!(tape::SparseVAFTape, op::SparseAffineOperation)
    pushfirst!(tape.operations, op)
    tape
end

# these aren't right, since they're mutating...
function MOIU.operate(::typeof(-), ::Type{T},
    tape::SparseVAFTape) where {T}
    d = MOI.output_dimension(tape)
    return add_operation!(tape, SparseAffineOperation(-sparse(1.0I, d, d), Zero(d)))
end

function MOIU.operate(::typeof(+), ::Type{T}, v::AbstractVector,
                            tape::SparseVAFTape) where {T}
    d = length(v)
    return add_operation!(tape, SparseAffineOperation(sparse(1.0I, d, d), v))
end

function MOIU.operate(::typeof(-), ::Type{T}, tape::SparseVAFTape,
    v::AbstractVector) where {T}
    d = length(v)
    return add_operation!(tape, SparseAffineOperation(sparse(1.0I, d, d), -v))
end

function MOIU.operate(::typeof(-), ::Type{T},
                    v::AbstractVector, tape::SparseVAFTape) where {T}
    d = length(v)
    return add_operation!(tape, SparseAffineOperation(sparse(-1.0I, d, d), v))
end

function MOIU.operate(::typeof(*), ::Type{T}, A::AbstractMatrix,
    tape::SparseVAFTape) where {T}
    return add_operation!(tape, SparseAffineOperation(A, Zero(size(A,1))))
end
####

USE_SPARSE() = true

function MOIU.operate(::typeof(*), ::Type{T}, A::AbstractMatrix,
                      v::MOI.VectorOfVariables) where {T}
    @assert size(A, 2) == length(v.variables)
    if USE_SPARSE()
        return SparseVAFTape([SparseAffineOperation(A, Zero(size(A,1)))], v.variables)

    else
        return VAFTape(tuple(AffineOperation(A, Zero(size(A,1)))), v.variables)
    end
end

function MOIU.operate(::typeof(+), ::Type{T}, v1::MOI.VectorOfVariables,
                      v2::MOI.ScalarAffineFunction) where {T}
    return MOIU.operate(+, T, scalar_fn(v1), v2)
end

function MOIU.operate(::typeof(+), ::Type{T}, v1::MOI.VectorAffineFunction,
                      v2::MOI.ScalarAffineFunction) where {T}
    return MOIU.operate(+, T, scalar_fn(v1), v2)
end

function MOIU.operate!(::typeof(+), ::Type{T}, v1::MOI.VectorAffineFunction,
                       v2::MOI.ScalarAffineFunction) where {T}
    return MOIU.operate!(+, T, scalar_fn(v1), v2)
end



# Other `MOIU.operate` methods

function MOIU.operate(::typeof(sum), ::Type{T}, vaf::MOI.VectorAffineFunction) where {T}
    return MOI.ScalarAffineFunction([vat.scalar_term for vat in vaf.terms],
                                    sum(vaf.constants))
end

function MOIU.operate(::typeof(sum), ::Type{T}, v::MOI.VectorOfVariables) where {T}
    return MOI.ScalarAffineFunction([MOI.ScalarAffineTerm{T}(one(T), vi)
                                     for vi in v.variables], zero(T))
end
