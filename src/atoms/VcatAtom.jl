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

function Base.getindex(
    x::VcatAtom,
    rows::AbstractVector{<:Real},
    cols::AbstractVector{<:Real},
)
    idx = 0
    rows = collect(rows) # make a mutable copy
    keep_children = ()
    for c in x.children
        # here are the row indices into `x` that point to `c`
        I = idx .+ (1:size(c, 1))
        if issubset(rows, I)
            # if all the row indices we want are in this one child, we can early exit
            if rows == I && cols == 1:size(c, 2)
                return c
            else
                return c[rows.-idx, cols]
            end
        elseif !isdisjoint(rows, I)
            # we have some but not all rows in this child, so keep it
            keep_children = (keep_children..., c)
            idx += size(c, 1)
        else # we can drop this child!
            # let's update `rows` to account for the removal
            l = last(I)
            for i in eachindex(rows)
                if rows[i] >= l
                    rows[i] -= length(I)
                end
            end
        end
    end
    # If we are here, the indices span multiple children.
    # We can't necessarily index each separately, since they may be out of order.
    # So we will defer to an `IndexAtom` on the remaining children
    remaining = VcatAtom(keep_children...)
    return IndexAtom(remaining, rows, cols)
end

function Base.getindex(x::VcatAtom, rows::AbstractVector{<:Real}, ::Colon)
    return x[rows, 1:size(x, 2)]
end

function Base.getindex(x::VcatAtom, ::Colon, cols::AbstractVector{<:Real})
    return x[1:size(x, 1), cols]
end

# linear indexing
# very similar to row-indexing above, but with linear indices
function Base.getindex(x::VcatAtom, inds::AbstractVector{<:Real})
    idx = 0
    inds = collect(inds)
    keep_children = ()
    for c in x.children
        I = idx .+ (1:length(c))
        if issubset(inds, I)
            if inds == I
                return c
            else
                return c[inds.-idx]
            end
        elseif !isdisjoint(inds, I)
            keep_children = (keep_children..., c)
            idx += length(c)
        else
            l = last(I)
            for i in eachindex(inds)
                if inds[i] >= l
                    inds[i] -= length(I)
                end
            end
        end
    end
    remaining = VcatAtom(keep_children...)
    return IndexAtom(remaining, inds)
end
