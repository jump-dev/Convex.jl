#############################################################################
# norm2.jl
# Handles the euclidean norm (also called frobenius norm for matrices)
# All expressions and atoms are subtpyes of AbstractExpr.
# Please read expressions.jl first.
#############################################################################
import LinearAlgebra.norm2
export EucNormAtom, norm2
export sign, monotonicity, curvature, conic_form!


struct EucNormAtom <: AbstractExpr
    head::Symbol
    id_hash::UInt64
    children::Tuple{AbstractExpr}
    size::Tuple{Int, Int}

    function EucNormAtom(x::AbstractExpr)
        children = (x,)
        return new(:norm2, hash(children), children, (1, 1))
    end
end

function sign(x::EucNormAtom)
    return Positive()
end

function monotonicity(x::EucNormAtom)
    return (sign(x.children[1]) * Nondecreasing(),)
end

function curvature(x::EucNormAtom)
    return ConvexVexity()
end

function evaluate(x::EucNormAtom)
    return norm(vec(evaluate(x.children[1])))
end


## Create a new variable euc_norm to represent the norm
## Additionally, create the second order conic constraint (euc_norm, x) in SOC
function conic_form!(x::EucNormAtom, unique_conic_forms::UniqueConicForms=UniqueConicForms())
    if !has_conic_form(unique_conic_forms, x)
        euc_norm = Variable()
        objective = conic_form!(euc_norm, unique_conic_forms)
        conic_form!(SOCConstraint(euc_norm, x.children[1]), unique_conic_forms)
        cache_conic_form!(unique_conic_forms, x, objective)
    end
    return get_conic_form(unique_conic_forms, x)
end

function norm2(x::AbstractExpr)
    if sign(x) == ComplexSign()
        return EucNormAtom([real(x);imag(x)])
    else
        return EucNormAtom(x)
    end
end
