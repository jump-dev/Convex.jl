import Base.sqrt
export ScaledGeoMeanAtom, scaledgeomean, sqrt
export sign, monotonicity, curvature, conic_form!

function power_of_2_gt(n::Int)
    Int(2^ceil(log(2,n)))
end

struct ScaledGeoMeanAtom <: AbstractExpr
    head::Symbol
    id_hash::UInt64
    children::Tuple{AbstractExpr, AbstractExpr}
    size::Tuple{Int, Int}

    function ScaledGeoMeanAtom(x::AbstractExpr)
        children = (x)
        return new(:scaledgeomean, hash(children), children, x.size)
    end
end

function sign(q::ScaledGeoMeanAtom)
    return Positive()
end

function monotonicity(q::ScaledGeoMeanAtom)
    return (Nondecreasing(),)
end

function curvature(q::ScaledGeoMeanAtom)
    return ConcaveVexity()
end

function evaluate(q::ScaledGeoMeanAtom)
    nbar = power_of_2_gt(length(q.children))
    return prod(evaluate(q.children[1]))^(1/nbar)
end

function conic_form!(q::ScaledGeoMeanAtom, unique_conic_forms::UniqueConicForms=UniqueConicForms())
    if !has_conic_form(unique_conic_forms, q)
        t = Variable()
        qol_objective = conic_form!(t, unique_conic_forms)
        x, y = q.children
        conic_form!(SOCElemConstraint(y + x, y - x, 2 * t), unique_conic_forms)
        conic_form!(x >= 0, unique_conic_forms)
        conic_form!(y >= 0, unique_conic_forms)
        cache_conic_form!(unique_conic_forms, q, qol_objective)
    end
    return get_conic_form(unique_conic_forms, q)
end

function scaledgeomean(x::AbstractExpr)
    if length(x) > 2
        return ScaledGeoMeanAtom(x)
    elseif length(x) == 2
        return geomean(x[1],x[2])
    else
        return x
    end
end
