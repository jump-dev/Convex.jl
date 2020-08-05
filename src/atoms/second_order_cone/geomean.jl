import Base.sqrt

struct GeoMeanAtom <: AbstractExpr
    head::Symbol
    id_hash::UInt64
    children::NTuple{N, AbstractExpr} where N
    size::Tuple{Int, Int}

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
    return prod.(evaluate.(q.children)).^Ref(1/n)
end

function template(q::GeoMeanAtom, context::Context{T}) where T
    n = length(q.children)
    x = first(q.children)
    t = Variable(size(x))
    for i = 1:size(x,1), j = 1:size(x,2)
        f = operate(vcat, T, template(t[i,j], context), ( template(y[i,j], context) for y in q.children)...)
        MOI_add_constraint(context.model, f, MOI.GeometricMeanCone(n+1))
    end
    return template(t, context)
end

geomean(x::AbstractExpr, y::AbstractExpr) = GeoMeanAtom(x, y)
sqrt(x::AbstractExpr) = GeoMeanAtom(x, Constant(ones(x.size[1], x.size[2])))
