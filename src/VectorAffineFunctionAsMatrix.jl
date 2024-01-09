# This is a variant of MOI.VectorAffineFunction which represents
# the transformation `matrix * variables + vector` lazily.
struct VectorAffineFunctionAsMatrix{T,V}
    aff::SparseAffineOperation{T}
    variables::V
end

function Base.isequal(
    A::VectorAffineFunctionAsMatrix,
    B::VectorAffineFunctionAsMatrix,
)
    return isequal(A.variables, B.variables) && isequal(A.aff, B.aff)
end

function MOI.output_dimension(v::VectorAffineFunctionAsMatrix)
    return size(v.matrix, 1)
end

function to_vaf(tape::SparseTape)
    op = SparseAffineOperation(tape)
    return to_vaf(VectorAffineFunctionAsMatrix(op, tape.variables))
end

# convert to a usual VAF
function to_vaf(vaf_as_matrix::VectorAffineFunctionAsMatrix{T}) where {T}
    I, J, V = SparseArrays.findnz(vaf_as_matrix.aff.matrix)
    vats = Vector{MOI.VectorAffineTerm{T}}(undef, length(I))
    for (idx, n) in enumerate(eachindex(I, J, V))
        i = I[n]
        j = J[n]
        v = V[n]

        vats[idx] = MOI.VectorAffineTerm{T}(
            i,
            MOI.ScalarAffineTerm{T}(v, vaf_as_matrix.variables[j]),
        )
    end

    return MOI.VectorAffineFunction{T}(vats, vaf_as_matrix.aff.vector)
end

# method for adding constraints and coverting to standard VAFs as needed
function MOI_add_constraint(model, f, set)
    return MOI.add_constraint(model, f, set)
end

function MOI_add_constraint(
    model,
    f::Union{VectorAffineFunctionAsMatrix,SparseTape},
    set,
)
    return MOI_add_constraint(model, to_vaf(f), set)
end
