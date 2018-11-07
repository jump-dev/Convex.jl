import Base.conj
export conj
export sign, curvature, monotonicity, evaluate, conic_form!
struct ConjugateAtom <: AbstractExpr
    head::Symbol
    id_hash::UInt64
    children::Tuple{AbstractExpr}
    size::Tuple{Int, Int}

    function ConjugateAtom(x::AbstractExpr)
        children = (x,)
        return new(:conj, hash(children), children, (x.size[1], x.size[2]))
    end
end

function sign(x::ConjugateAtom)
    return sign(x.children[1])
end

function monotonicity(x::ConjugateAtom)
    return (Nondecreasing(),)
end

function curvature(x::ConjugateAtom)
    return ConstVexity()
end

function evaluate(x::ConjugateAtom)
    return conj(evaluate(x.children[1]))
end

function conic_form!(x::ConjugateAtom, unique_conic_forms::UniqueConicForms=UniqueConicForms())
    if !has_conic_form(unique_conic_forms, x)
        objective = conic_form!(x.children[1], unique_conic_forms)
        for var in keys(objective)
            x1 = conj(objective[var][1])
            x2 = conj(objective[var][2])
            objective[var] = (x1,x2)
        end
        cache_conic_form!(unique_conic_forms, x, objective)
    end
    return get_conic_form(unique_conic_forms, x)
end


conj(x::AbstractExpr) = ConjugateAtom(x)
conj(x::Constant) = Constant(conj(x.value))
