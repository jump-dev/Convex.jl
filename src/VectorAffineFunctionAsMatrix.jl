# This is a variant of MOI.VectorAffineFunction which represents
# the transformation `matrix * variables + vector` lazily.
struct VectorAffineFunctionAsMatrix{M,B,V}
    matrix::M
    vector::B
    variables::V
end

struct AffineOperation{M, V}
    matrix::M
    vector::V
end

struct VAFTape{T <: Tuple}
    operations::T
    variables::Vector{MOI.VariableIndex}
end

# A simple type representing a vector of zeros. Maybe should include the size or use FillArrays or similar.
struct Zero end

Base.:(+)(a, ::Zero) = a
Base.:(+)(::Zero, a) = a
Base.:(-)(::Zero, a) = -a
Base.:(-)(a, ::Zero) = a
Base.:(-)(::Zero) = Zero()
Base.:(*)(A, ::Zero) = Zero()

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

function MOI.output_dimension(v::VectorAffineFunctionAsMatrix)
    return size(v.matrix, 1)
end

function popfirst(t::Tuple)
    return first(t), Base.tail(t)
end

# the order here is chosen deliberately: generally, the final output is 1-dimensional,
# so we multiply from left to right.
function compile!(tape::VAFTape)# -> AffineOperation
    t = tape.operations
    op, t = popfirst(t)
    # we will modify op in place
    mat = op.matrix
    vec = op.vector
    while !isempty(t)
        newop, t = popfirst(t)
        # op.vector += op.matrix * newop.vector
        vec = try_mul!!(vec, op.matrix, newop.vector, 1, 1)

        # mat = mat * newop.matrix
        mat = try_mul!!(mat, newop.matrix)
    end
    return AffineOperation(mat, vec)
end

"""
    try_mul!!(A, B)

Performs `A*B` and returns the result, possibly overwriting either `A` or `B`.
"""
function try_mul!! end

function try_mul!!(C, A, B, α, β)
    mul!(C, A, B, α, β)
end

function try_mul!!(C, A, ::Zero, α, β)
    try_mul!!(C, β)
end

function try_mul!!(A, B)
    rmul!(A, B)
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

function to_vaf(tape::VAFTape)
    op = compile!(tape)
    to_vaf(VectorAffineFunctionAsMatrix(op.matrix, op.vector, tape.variables))
end


# convert to a usual VAF
function to_vaf(vaf_as_matrix::VectorAffineFunctionAsMatrix{<:SparseMatrixCSC}, context)
    T = eltype(vaf_as_matrix.matrix)
    I, J, V = findnz(vaf_as_matrix.matrix)
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

    return MOI.VectorAffineFunction{T}(vats, vaf_as_matrix.vector)
end

function to_vaf(vaf_as_matrix::VectorAffineFunctionAsMatrix{<:AbstractMatrix})
    T = eltype(vaf_as_matrix.matrix)
    vats = MOI.VectorAffineTerm{T}[]
    M = vaf_as_matrix.matrix
    for i in 1:size(M, 1)
        for j in 1:size(M, 2)
            # this check turns out to be key for performance on the test problem
            # which means I suspect something is being densely represented when possibly it should be sparse
            iszero(M[i, j]) && continue
            push!(vats,
                  MOI.VectorAffineTerm{T}(i,
                                          MOI.ScalarAffineTerm{T}(M[i, j],
                                                                  vaf_as_matrix.variables[j])))
        end
    end
    return MOI.VectorAffineFunction{T}(vats, vaf_as_matrix.vector)
end

# method for adding constraints and coverting to standard VAFs as needed
MOI_add_constraint(model, f, set) = MOI.add_constraint(model, f, set)

function MOI_add_constraint(model, f::VectorAffineFunctionAsMatrix, set)
    return MOI.add_constraint(model, to_vaf(f), set)
end

function MOI_add_constraint(model, f::VAFTape, set)
    return MOI.add_constraint(model, to_vaf(f), set)
end
# # `MOIU.operate` methods for `VAFTape`
# function MOIU.operate(f::F, ::Type{T}, tape::VAFTape, args) where {F, T}
#     return VAFTape((VAFOperation(f, T, args), tape.operations...), tape.vafasmatrix)
# end

add_operation(op::AffineOperation, tape::VAFTape) = VAFTape((op, tape.operations...), tape.variables)

function MOIU.operate(::typeof(-), ::Type{T},
    tape::VAFTape) where {T}
    return add_operation(AffineOperation(-I, Zero()), tape)
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
                      v::MOI.VectorOfVariables) where {T}
    @assert size(A, 2) == length(v.variables)
    return VAFTape(tuple(AffineOperation(A, Zero())), v.variables)
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

function MOIU.operate(::typeof(*), ::Type{T}, A::AbstractMatrix,
                      tape::VAFTape) where {T}
    return add_operation(AffineOperation(A, Zero()), tape)
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
