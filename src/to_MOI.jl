# Copyright (c) 2014: Madeleine Udell and contributors
#
# Use of this source code is governed by a BSD-style license that can be found
# in the LICENSE file or at https://opensource.org/license/bsd-2-clause

# This is a variant of MOI.VectorAffineFunction which represents
# the transformation `matrix * variables + vector` lazily.
to_vaf(tape::SPARSE_VECTOR) = tape

# convert to a usual VAF
function to_vaf(tape::SparseTape{T}) where {T}
    I, J, V = SparseArrays.findnz(tape.operation.matrix)
    vats = Vector{MOI.VectorAffineTerm{T}}(undef, length(I))
    for (idx, n) in enumerate(eachindex(I, J, V))
        i = I[n]
        j = J[n]
        v = V[n]

        vats[idx] = MOI.VectorAffineTerm{T}(
            i,
            MOI.ScalarAffineTerm{T}(v, tape.variables[j]),
        )
    end

    return MOI.VectorAffineFunction{T}(vats, tape.operation.vector)
end

# method for adding constraints and converting to standard VAFs as needed
function MOI_add_constraint(model, f, set)
    return MOI.add_constraint(model, f, set)
end

function MOI_add_constraint(model, f::SPARSE_VECTOR{T}, set) where {T}
    g = MOI.VectorAffineFunction{T}(MOI.VectorAffineTerm{T}[], f)
    return MOI.add_constraint(model, g, set)
end

function MOI_add_constraint(model, f::SparseTape, set)
    return MOI_add_constraint(model, to_vaf(f), set)
end
