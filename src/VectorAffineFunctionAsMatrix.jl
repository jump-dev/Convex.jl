# Copyright (c) 2014: Madeleine Udell and contributors
#
# Use of this source code is governed by a BSD-style license that can be found
# in the LICENSE file or at https://opensource.org/license/bsd-2-clause

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

function to_saf(tape::SparseTape{T}, output_index) where {T}
    A = tape.operation.matrix
    m, n = size(A)
    sats = MOI.ScalarAffineTerm{T}[]
    rows = SparseArrays.rowvals(A)
    vals = SparseArrays.nonzeros(A)
    for j = 1:n # for each column
        for i in SparseArrays.nzrange(A, j)
            row = rows[i]
            row == output_index || continue
            val = vals[i]
            push!(sats, MOI.ScalarAffineTerm{T}(val, tape.variables[j]))
        end
    end
    return MOI.ScalarAffineFunction{T}(sats, tape.operation.vector[output_index])
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
