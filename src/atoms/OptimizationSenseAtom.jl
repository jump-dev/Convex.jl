# Copyright (c) 2014: Madeleine Udell and contributors
#
# Use of this source code is governed by a BSD-style license that can be found
# in the LICENSE file or at https://opensource.org/license/bsd-2-clause

mutable struct OptimizationSenseAtom <: AbstractExpr
    children::Tuple{AbstractExpr}
    sense::MOI.OptimizationSense
    size::Tuple{Int,Int}

    function OptimizationSenseAtom(
        x::AbstractExpr,
        sense::MOI.OptimizationSense,
    )
        @assert sense != MOI.FEASIBILITY_SENSE
        return new((x,), sense, size(x))
    end
end

head(io::IO, x::OptimizationSenseAtom) = print(io, x.sense)

Base.sign(x::OptimizationSenseAtom) = sign(only(x.children))

monotonicity(x::OptimizationSenseAtom) = (Nondecreasing(),)

function curvature(x::OptimizationSenseAtom)
    if x.sense == MOI.MIN_SENSE
        return ConvexVexity()
    else
        @assert x.sense == MOI.MAX_SENSE
        return ConcaveVexity()
    end
end

evaluate(x::OptimizationSenseAtom) = evaluate(only(x.children))

function new_conic_form!(context::Context, x::OptimizationSenseAtom)
    return new_conic_form!(context, only(x.children))
end
