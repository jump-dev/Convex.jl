import Base.sqrt

mutable struct GeoMeanAtom <: AbstractExpr
    children::NTuple{N,AbstractExpr} where {N}
    size::Tuple{Int,Int}

    function GeoMeanAtom(args::AbstractExprOrValue...)
        args = Tuple(arg isa Value ? constant(arg) : arg for arg in args)
        sz = size(first(args))
        if any(!=(sz), size.(args))
            error("geo mean must take arguments of the same size")
        elseif any(x -> sign(x) == ComplexSign(), args)
            error("The arguments should be real, not complex")
        else
            children = args
            return new(children, sz)
        end
    end
end

head(io::IO, ::GeoMeanAtom) = print(io, "geomean")

function sign(q::GeoMeanAtom)
    return Positive()
end

function monotonicity(q::GeoMeanAtom)
    return fill(Nondecreasing(), length(q.children))
end

function curvature(q::GeoMeanAtom)
    return ConcaveVexity()
end

_geomean(scalar_args...) = prod(scalar_args)^(1 / length(scalar_args))
function evaluate(q::GeoMeanAtom)
    n = length(q.children)
    children = evaluate.(q.children)
    c1 = first(children)
    return [_geomean((children[i][I] for i in 1:n)...) for I in eachindex(c1)]
end

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
sqrt(x::AbstractExpr) = GeoMeanAtom(x, constant(ones(x.size[1], x.size[2])))
