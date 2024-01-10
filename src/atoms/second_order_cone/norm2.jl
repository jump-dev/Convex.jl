mutable struct EucNormAtom <: AbstractExpr
    children::Tuple{AbstractExpr}
    size::Tuple{Int,Int}

    EucNormAtom(x::AbstractExpr) = new((x,), (1, 1))
end

head(io::IO, ::EucNormAtom) = print(io, "norm2")

Base.sign(::EucNormAtom) = Positive()

monotonicity(x::EucNormAtom) = (sign(x.children[1]) * Nondecreasing(),)

curvature(::EucNormAtom) = ConvexVexity()

evaluate(x::EucNormAtom) = norm(vec(evaluate(x.children[1])))

function new_conic_form!(context::Context{T}, A::EucNormAtom) where {T}
    x = only(AbstractTrees.children(A))
    t = conic_form!(context, Variable())
    obj = conic_form!(context, x)
    f = operate(vcat, T, sign(A), t, obj)
    MOI_add_constraint(context.model, f, MOI.SecondOrderCone(length(x) + 1))
    return t
end

function LinearAlgebra.norm2(x::AbstractExpr)
    if sign(x) == ComplexSign()
        return EucNormAtom([real(x); imag(x)])
    end
    return EucNormAtom(x)
end
