# Copyright (c) 2014: Madeleine Udell and contributors
#
# Use of this source code is governed by a BSD-style license that can be found
# in the LICENSE file or at https://opensource.org/license/bsd-2-clause

_to_expr(x::AbstractExpr) = x
_to_expr(x::Value) = constant(x)

mutable struct GeoMeanAtom <: AbstractExpr
    children::Vector{AbstractExpr}
    size::Tuple{Int,Int}

    function GeoMeanAtom(args::Union{AbstractExpr,Value}...)
        new_args = AbstractExpr[_to_expr(arg) for arg in args]
        sz = size(first(new_args))
        if any(!=(sz), size.(new_args))
            error("[GeoMeanAtom] geomean must take arguments of the same size")
        elseif any(x -> sign(x) == ComplexSign(), new_args)
            error("[GeoMeanAtom] the arguments must be real, not complex")
        end
        return new(new_args, sz)
    end
end

head(io::IO, ::GeoMeanAtom) = print(io, "geomean")

Base.sign(::GeoMeanAtom) = Positive()

function monotonicity(q::GeoMeanAtom)
    return ntuple(_ -> Nondecreasing(), length(q.children))
end

curvature(::GeoMeanAtom) = ConcaveVexity()

_geomean(scalar_args...) = prod(scalar_args)^(1 / length(scalar_args))

evaluate(q::GeoMeanAtom) = _geomean.(evaluate.(q.children)...)

function new_conic_form!(context::Context{T}, q::GeoMeanAtom) where {T}
    n = length(q.children)
    x = first(q.children)
    t = Variable(size(x))
    t_tape = conic_form!(context, t)
    child_fns = map(q.children) do y
        return MOI.Utilities.scalarize(to_vaf(conic_form!(context, y)))
    end
    for (i, ti) in enumerate(t_tape.variables)
        f = MOI.Utilities.operate(vcat, T, ti, getindex.(child_fns, i)...)
        MOI.add_constraint(context.model, f, MOI.GeometricMeanCone(n + 1))
    end
    return t_tape
end

geomean(args::Union{AbstractExpr,Value}...) = GeoMeanAtom(args...)

function Base.sqrt(x::AbstractExpr)
    return GeoMeanAtom(x, constant(ones(x.size[1], x.size[2])))
end
