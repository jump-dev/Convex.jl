# This is a variant of MOI.VectorAffineFunction which represents
# the transformation `matrix * variables + vector` lazily.
struct VectorAffineFunctionAsMatrix{M,B,V}
    matrix::M
    vector::B
    variables::V
end

function MOI.output_dimension(v::VectorAffineFunctionAsMatrix)
    return size(v.matrix, 1)
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

# `MOIU.operate` methods for `VectorAffineFunctionAsMatrix`

function MOIU.operate(::typeof(-), ::Type{T},
                      vafasmatrix::VectorAffineFunctionAsMatrix) where {T}
    return VectorAffineFunctionAsMatrix(-vafasmatrix.matrix, -vafasmatrix.vector,
                                        vafasmatrix.variables)
end

function MOIU.operate(::typeof(+), ::Type{T}, v::AbstractVector,
                      vafasmatrix::VectorAffineFunctionAsMatrix) where {T}
    return VectorAffineFunctionAsMatrix(vafasmatrix.matrix, vafasmatrix.vector + v,
                                        vafasmatrix.variables)
end

function MOIU.operate(::typeof(-), ::Type{T}, vafasmatrix::VectorAffineFunctionAsMatrix,
                      v::AbstractVector) where {T}
    return VectorAffineFunctionAsMatrix(vafasmatrix.matrix, vafasmatrix.vector - v,
                                        vafasmatrix.variables)
end

# A simple type representing a vector of zeros. Maybe should include the size or use FillArrays or similar.
struct Zero end

Base.:(+)(a, z::Zero) = a
Base.:(+)(z::Zero, a) = a
Base.:(-)(z::Zero, a) = -a
Base.:(-)(a, z::Zero) = a
Base.:(-)(z::Zero) = z
Base.:(*)(A, z::Zero) = z

function MOIU.operate(::typeof(*), ::Type{T}, A::AbstractMatrix,
                      v::MOI.VectorOfVariables) where {T}
    @assert size(A, 2) == length(v.variables)
    return VectorAffineFunctionAsMatrix(A, Zero(), v.variables)
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
                      vafasmatrix::VectorAffineFunctionAsMatrix) where {T}
    return VectorAffineFunctionAsMatrix(A * vafasmatrix.matrix, A * vafasmatrix.vector,
                                        vafasmatrix.variables)
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
