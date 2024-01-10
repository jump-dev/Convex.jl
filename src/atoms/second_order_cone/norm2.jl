#############################################################################
# LinearAlgebra.norm2.jl
# Handles the euclidean norm (also called frobenius norm for matrices)
# All expressions and atoms are subtpyes of AbstractExpr.
# Please read expressions.jl first.
#############################################################################

mutable struct EucNormAtom <: AbstractExpr
    children::Tuple{AbstractExpr}
    size::Tuple{Int,Int}

    function EucNormAtom(x::AbstractExpr)
        children = (x,)
        return new(children, (1, 1))
    end
end

head(io::IO, ::EucNormAtom) = print(io, "norm2")

function Base.sign(x::EucNormAtom)
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
function new_conic_form!(context::Context{T}, A::EucNormAtom) where {T}
    obj = conic_form!(context, only(AbstractTrees.children(A)))

    x = only(AbstractTrees.children(A))
    d = length(x)

    t_obj = conic_form!(context, Variable())

    f = operate(vcat, T, sign(A), t_obj, obj)
    set = MOI.SecondOrderCone(d + 1)
    MOI_add_constraint(context.model, f, set)

    return t_obj
end

function LinearAlgebra.norm2(x::AbstractExpr)
    if sign(x) == ComplexSign()
        return EucNormAtom([real(x); imag(x)])
    else
        return EucNormAtom(x)
    end
end
