# Copyright (c) 2014: Madeleine Udell and contributors
#
# Use of this source code is governed by a BSD-style license that can be found
# in the LICENSE file or at https://opensource.org/license/bsd-2-clause

mutable struct ComplexTape{T}
    real_tape::SparseTape{T}
    imag_tape::SparseTape{T}

    function ComplexTape(re::SparseTape{T}, im::SparseTape{T}) where {T}
        MOI.output_dimension(re) == MOI.output_dimension(im) ||
            DimensionMismatch()
        return new{T}(re, im)
    end
end

MOI.output_dimension(c::ComplexTape) = MOI.output_dimension(c.real_tape) # same value as for imag
Base.real(c::ComplexTape) = c.real_tape
Base.imag(c::ComplexTape) = c.imag_tape

function MOI_add_constraint(
    model,
    f::ComplexTape,
    set::Union{MOI.Zeros,MOI.Nonnegatives,MOI.Nonpositives},
)
    re_inds = MOI_add_constraint(model, to_vaf(real(f)), set)
    im_inds = MOI_add_constraint(model, to_vaf(imag(f)), set)
    return (re_inds, im_inds)
end
