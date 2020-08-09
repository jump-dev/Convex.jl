# A simple type representing a vector of zeros. Maybe should include the size or use FillArrays or similar.
struct Zero
    len::Int
 end

Base.:(+)(a, ::Zero) = a
Base.:(+)(::Zero, a) = a
Base.:(+)(z1::Zero, z2::Zero) = (@assert z1.len == z2.len; z1)
Base.:(-)(::Zero, a) = -a
Base.:(-)(a, ::Zero) = a
Base.:(-)(z::Zero) = z
Base.:(*)(A, z::Zero) = Zero(size(A, 1))
Base.:(*)(::LinearAlgebra.UniformScaling, z::Zero) = z
Base.size(z::Zero) = (z.len,)
Base.length(z::Zero) = z.len
SparseArrays.sparse(z::Zero) = spzeros(z.len)
Base.convert(::Type{Vector{T}}, z::Zero) where {T} = zeros(T, z.len)

include("VAFTape.jl")
include("SparseVAFTape.jl")
include("SparseVAFTape2.jl")
include("SparseVAFTape3.jl")


const VAFTapes = Union{VAFTape, SparseVAFTape, SparseVAFTape2, SparseVAFTape3}

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
    op = AffineOperation(tape)
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

## this may be wrong
# function to_vaf_as_matrix(vaf::MOI.VectorAffineFunction{T}) where T
#     I = Int[]
#     J = Int[]
#     V = T[]
#     variables = MOI.VariableIndex[]
#     var_dict = Dict{MOI.VariableIndex, Int}()
#     for term in vaf.terms
#         push!(I, term.output_index)
#         push!(V, term.scalar_term.coefficient)
#         var = term.scalar_term.variable_index
#         if haskey(var_dict, var)
#             j = var_dict[var]
#         else
#             push!(variables, var)
#             var_dict[var] = j = length(variables)
#         end
#         push!(J, j)
#     end
#     n = length(variables)
#     m = MOI.output_dimension(vaf)
#     return VectorAffineFunctionAsMatrix(AffineOperation(sparse(I, J, V, m, n), vaf.constants), variables)
# end


function to_vaf(vaf_as_matrix::VectorAffineFunctionAsMatrix{<:UniformScaling})
    T = eltype(vaf_as_matrix.aff.matrix)
    vats = MOI.VectorAffineTerm{T}[]
    x = I.λ
    for i in 1:length(vaf_as_matrix.variables)
        push!(vats,
                MOI.VectorAffineTerm{T}(i,
                                        MOI.ScalarAffineTerm{T}(x,
                                                                vaf_as_matrix.variables[i])))
    end
    return MOI.VectorAffineFunction{T}(vats, vaf_as_matrix.aff.vector)
end

# method for adding constraints and coverting to standard VAFs as needed
function MOI_add_constraint(model, f, set)
    return MOI.add_constraint(model, f, set)
end

function MOI_add_constraint(model, f::Union{VectorAffineFunctionAsMatrix, VAFTapes}, set)
    return MOI_add_constraint(model, to_vaf(f), set)
end
