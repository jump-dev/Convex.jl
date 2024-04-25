# Copyright (c) 2014: Madeleine Udell and contributors
#
# Use of this source code is governed by a BSD-style license that can be found
# in the LICENSE file or at https://opensource.org/license/bsd-2-clause

mutable struct HcatAtom <: AbstractExpr
    children::Tuple
    size::Tuple{Int,Int}

    function HcatAtom(args...)
        args = convert.(AbstractExpr, args)
        num_cols, num_rows = 0, args[1].size[1]
        for arg in args
            if arg.size[1] != num_rows
                msg = "[HcatAtom] cannot stack expressions of incompatible size. Got $(arg.size[1]) expected $num_rows."
                throw(DimensionMismatch(msg))
            end
            num_cols += arg.size[2]
        end
        return new(args, (num_rows, num_cols))
    end
end

head(io::IO, ::HcatAtom) = print(io, "hcat")

Base.sign(x::HcatAtom) = sum(map(sign, x.children))

monotonicity(x::HcatAtom) = ntuple(_ -> Nondecreasing(), length(x.children))

curvature(::HcatAtom) = ConstVexity()

evaluate(x::HcatAtom) = hcat(map(evaluate, x.children)...)

function new_conic_form!(context::Context{T}, x::HcatAtom) where {T}
    objectives = map(c -> conic_form!(context, c), AbstractTrees.children(x))
    # Suppose the child objectives for two children e1 (2 x 1) and e2 (2 x 2)
    # look something like
    #  e1: x => 1 2 3
    #           4 5 6
    #      y => 2 4
    #           7 8
    #  e2: x => 1 1 1
    #           2 2 2
    #           3 3 3
    #           4 4 4
    # The objective of [e1 e2] will look like
    #            x => 1 2 3
    #                 4 5 6
    #                 1 1 1
    #                 2 2 2
    #                 3 3 3
    #                 4 4 4
    #            y => 2 4
    #                 7 8
    #                 0 0
    #                 0 0
    #                 0 0
    #                 0 0
    # builds the objective by aggregating a list of coefficients for each
    # variable from each child objective, and then vertically concatenating them
    return operate(vcat, T, sign(x), objectives...)
end
# TODO: fix piracy!

# * `Value` is not owned by Convex.jl
# * splatting creates zero-argument functions, which again are not owned by Convex.jl
Base.hcat(args::AbstractExpr...) = HcatAtom(args...)

function Base.hcat(args::Union{AbstractExpr,Value}...)
    if all(Base.Fix2(isa, Value), args)
        return Base.cat(args..., dims = Val(2))
    end
    return HcatAtom(args...)
end

# TODO: implement vertical concatenation in a more efficient way
Base.vcat(args::AbstractExpr...) = transpose(HcatAtom(map(transpose, args)...))

function Base.vcat(args::Union{AbstractExpr,Value}...)
    if all(Base.Fix2(isa, Value), args)
        return Base.cat(args..., dims = Val(1))
    end
    return transpose(HcatAtom(map(transpose, args)...))
end

function Base.hvcat(
    rows::Tuple{Vararg{Int}},
    args::Union{AbstractExpr,Value}...,
)
    nbr = length(rows)
    rs = Vector{Any}(undef, nbr)
    a = 1
    for i in 1:nbr
        rs[i] = HcatAtom(args[a:a-1+rows[i]]...)
        a += rows[i]
    end
    return vcat(rs...)
end
