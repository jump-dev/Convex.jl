import Base.sqrt

struct GeoMeanAtom <: AbstractExpr
    head::Symbol
    id_hash::UInt64
    children::NTuple{N,AbstractExpr} where {N}
    size::Tuple{Int,Int}

    function GeoMeanAtom(args::AbstractExpr...)
        sz = size(first(args))
        if any(!=(sz), size.(args))
            error("geo mean must take arguments of the same size")
        elseif any(x -> sign(x) == ComplexSign(), args)
            error("The arguments should be real, not complex")
        else
            children = args
            return new(:geomean, hash(children), children, sz)
        end
    end
end

function sign(q::GeoMeanAtom)
    return Positive()
end

function monotonicity(q::GeoMeanAtom)
    return (Nondecreasing(), Nondecreasing())
end

function curvature(q::GeoMeanAtom)
    return ConcaveVexity()
end

function evaluate(q::GeoMeanAtom)
    n = length(q.children)
    return prod.(evaluate.(q.children)) .^ Ref(1 / n)
end

function conic_form!(context::Context{T}, q::GeoMeanAtom) where {T}
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

geomean(x::AbstractExpr, y::AbstractExpr) = GeoMeanAtom(x, y)
sqrt(x::AbstractExpr) = GeoMeanAtom(x, constant(ones(x.size[1], x.size[2])))
