_to_expr(x::AbstractExpr) = x
_to_expr(x::Value) = constant(x)

mutable struct GeoMeanAtom <: AbstractExpr
    children::Vector{AbstractExpr}
    size::Tuple{Int,Int}

    function GeoMeanAtom(args::AbstractExprOrValue...)
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
    for i in 1:size(x, 1), j in 1:size(x, 2)
        f = operate(
            vcat,
            T,
            sign(q),
            conic_form!(context, t[i, j]),
            (conic_form!(context, y[i, j]) for y in q.children)...,
        )
        MOI_add_constraint(context.model, f, MOI.GeometricMeanCone(n + 1))
    end
    return conic_form!(context, t)
end

geomean(args::AbstractExprOrValue...) = GeoMeanAtom(args...)

function Base.sqrt(x::AbstractExpr)
    return GeoMeanAtom(x, constant(ones(x.size[1], x.size[2])))
end
