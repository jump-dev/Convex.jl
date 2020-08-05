#############################################################################
# norm2.jl
# Handles the euclidean norm (also called frobenius norm for matrices)
# All expressions and atoms are subtpyes of AbstractExpr.
# Please read expressions.jl first.
#############################################################################
import LinearAlgebra.norm2


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

function template(A::EucNormAtom, context::Context{T}) where T
    obj = template(only(children(A)), context)

    x = only(children(A))
    d = length(x)

    t_obj = template(Variable(), context)

    f = operate(vcat, T, t_obj, obj)
    set = MOI.SecondOrderCone(d + 1)
    MOI_add_constraint(context.model, f, set)

    return t_obj
end

function norm2(x::AbstractExpr)
    if sign(x) == ComplexSign()
        return EucNormAtom([real(x);imag(x)])
    else
        return EucNormAtom(x)
    end
end
