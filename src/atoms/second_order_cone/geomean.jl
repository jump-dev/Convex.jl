import Base.sqrt
export GeoMeanAtom, geomean, sqrt
export sign, monotonicity, curvature, conic_form!

struct GeoMeanAtom <: AbstractExpr
    head::Symbol
    id_hash::UInt64
    children::Tuple{AbstractExpr, AbstractExpr}
    size::Tuple{Int, Int}

    function GeoMeanAtom(x::AbstractExpr, y::AbstractExpr)
        if x.size != y.size
            error("geo mean must take two arguments of the same size")
        elseif sign(x) == ComplexSign() || sign(y) == ComplexSign()
            error("Both the arguments should be real instead they are $(sign(x)) and $(sign(y))")
        else
            children = (x, y)
            return new(:geomean, hash(children), children, x.size)
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
    return sqrt.(evaluate(q.children[1]) .* evaluate(q.children[2]))
end

function conic_form!(q::GeoMeanAtom, unique_conic_forms::UniqueConicForms=UniqueConicForms())
    if !has_conic_form(unique_conic_forms, q)
        sz = q.children[1].size
        t = Variable(sz[1], sz[2])
        qol_objective = conic_form!(t, unique_conic_forms)
        x, y = q.children
        conic_form!(SOCElemConstraint(y + x, y - x, 2 * t), unique_conic_forms)
        conic_form!(x >= 0, unique_conic_forms)
        conic_form!(y >= 0, unique_conic_forms)
        cache_conic_form!(unique_conic_forms, q, qol_objective)
    end
    return get_conic_form(unique_conic_forms, q)
end

geomean(x::AbstractExpr, y::AbstractExpr) = GeoMeanAtom(x, y)
sqrt(x::AbstractExpr) = GeoMeanAtom(x, Constant(Ones(x.size[1], x.size[2])))
