# Copyright (c) 2014: Madeleine Udell and contributors
#
# Use of this source code is governed by a BSD-style license that can be found
# in the LICENSE file or at https://opensource.org/license/bsd-2-clause

mutable struct VcatAtom <: AbstractExpr
    children::Tuple
    size::Tuple{Int,Int}

    function VcatAtom(args...)
        args = convert.(AbstractExpr, args)
        num_rows, num_cols = 0, args[1].size[2]
        for arg in args
            if arg.size[2] != num_cols
                msg = "[VcatAtom] cannot stack expressions of incompatible size. Got $(arg.size[2]) expected $num_cols."
                throw(DimensionMismatch(msg))
            end
            num_rows += arg.size[1]
        end
        return new(args, (num_rows, num_cols))
    end
end

head(io::IO, ::VcatAtom) = print(io, "vcat")

Base.sign(x::VcatAtom) = sum(map(sign, x.children))

monotonicity(x::VcatAtom) = ntuple(_ -> Nondecreasing(), length(x.children))

curvature(::VcatAtom) = ConstVexity()

evaluate(x::VcatAtom) = reduce(vcat, collect(map(evaluate, x.children)))

function new_conic_form!(context::Context{T}, x::VcatAtom) where {T}
    # Converting a VcatAtom to conic form is non-trivial. Consider two matrices:
    #   x = [1 3; 2 4]
    #   y = [5 7; 6 8]
    # with VcatAtom(x, y). The desired outcome is [1, 2, 5, 6, 3, 4, 7, 8].
    # If we naively convert the children to conic form and then vcat, we will
    # get:
    #   vcat([1, 2, 3, 4], [5, 6, 7, 8]) = [1, 2, 3, 4, 5, 6, 7, 8]
    # which is not what we are after. We need to first transpose each child to
    # get:
    #   x^T, y^T = [1 2; 3 4], [5 6; 7 8])
    # then hcat them to get:
    #   hcat(x^T, y^T) = [1 2 5 6; 3 4 7 8]
    # then transpose this to get:
    #   hcat(x^T, y^T)^T = [1 3; 2 4; 5 7; 6 8]
    # so our final conic form produces the desired
    #   [1, 2, 5, 6, 3, 4, 7, 8]
    return conic_form!(context, transpose(reduce(hcat, transpose.(x.children))))
end

Base.vcat(args::AbstractExpr...) = VcatAtom(args...)

function Base.vcat(args::Union{AbstractExpr,Value}...)
    if all(Base.Fix2(isa, Value), args)
        return Base.cat(args..., dims = Val(1))
    end
    return VcatAtom(args...)
end
