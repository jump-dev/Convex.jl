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


## Create a new variable euc_norm to represent the norm
## Additionally, create the second order conic constraint (euc_norm, x) in SOC
function conic_form!(x::EucNormAtom, unique_conic_forms::UniqueConicForms)
    if !has_conic_form(unique_conic_forms, x)
        euc_norm = Variable()
        objective = conic_form!(euc_norm, unique_conic_forms)
        conic_form!(SOCConstraint(euc_norm, x.children[1]), unique_conic_forms)
        cache_conic_form!(unique_conic_forms, x, objective)
    end
    return get_conic_form(unique_conic_forms, x)
end


function template(A::EucNormAtom, context::Context{T}) where T
    obj = template(only(children(A)), context)

    x = only(children(A))
    d = length(x)

    t_obj = template(Variable(), context)

    t_single_var_obj = MOI.SingleVariable(only(t_obj.variables))

    # we're going to convert early since we haven't defined `vcat`...
    if obj isa VectorAffineFunctionAsMatrix || obj isa VAFTapes
        obj = to_vaf(obj)
    end

    f = MOIU.operate(vcat, T, t_single_var_obj, obj)
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
