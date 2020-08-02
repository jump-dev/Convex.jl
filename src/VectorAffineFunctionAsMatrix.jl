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

include("VAFTape.jl")
include("SparseVAFTape.jl")
const VAFTapes = Union{VAFTape, SparseVAFTape}

# This is a variant of MOI.VectorAffineFunction which represents
# the transformation `matrix * variables + vector` lazily.
struct VectorAffineFunctionAsMatrix{M,B,V}
    aff::AffineOperation{M, B}
    variables::V
end

function Base.isequal(A::VectorAffineFunctionAsMatrix, B::VectorAffineFunctionAsMatrix)
    isequal(A.variables, B.variables) && isequal(A.aff, B.aff)
end




function MOI.output_dimension(v::VectorAffineFunctionAsMatrix)
    return size(v.matrix, 1)
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


USE_SPARSE() = true
# This is the entry point into the `VAFTape`-land.
# `USE_SPARSE` governs which path is taken.
function MOIU.operate(::typeof(*), ::Type{T}, A::AbstractMatrix,
                      v::MOI.VectorOfVariables) where {T}
    @assert size(A, 2) == length(v.variables)
    if USE_SPARSE()
        return SparseVAFTape([SparseAffineOperation(A, Zero(size(A,1)))], v.variables)

    else
        return VAFTape(tuple(AffineOperation(A, Zero(size(A,1)))), v.variables)
    end
end


# Other `MOIU.operate` methods

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

function MOIU.operate(::typeof(sum), ::Type{T}, vaf::MOI.VectorAffineFunction) where {T}
    return MOI.ScalarAffineFunction([vat.scalar_term for vat in vaf.terms],
                                    sum(vaf.constants))
end

function MOIU.operate(::typeof(sum), ::Type{T}, v::MOI.VectorOfVariables) where {T}
    return MOI.ScalarAffineFunction([MOI.ScalarAffineTerm{T}(one(T), vi)
                                     for vi in v.variables], zero(T))
end
