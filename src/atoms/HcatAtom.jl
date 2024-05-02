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

evaluate(x::HcatAtom) = reduce(hcat, collect(map(evaluate, x.children)))

function new_conic_form!(context::Context{T}, x::HcatAtom) where {T}
    args = map(c -> conic_form!(context, c), AbstractTrees.children(x))
    # MOI represents matrices by concatenating their columns, so even though
    # this is an HcatAtom, we built the conic form by vcat'ing the arguments.
    return operate(vcat, T, sign(x), args...)
end
