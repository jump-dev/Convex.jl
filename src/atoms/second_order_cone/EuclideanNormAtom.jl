mutable struct EuclideanNormAtom <: AbstractExpr
    children::Tuple{AbstractExpr}
    size::Tuple{Int,Int}

    EuclideanNormAtom(x::AbstractExpr) = new((x,), (1, 1))
end

head(io::IO, ::EuclideanNormAtom) = print(io, "norm2")

Base.sign(::EuclideanNormAtom) = Positive()

monotonicity(x::EuclideanNormAtom) = (sign(x.children[1]) * Nondecreasing(),)

curvature(::EuclideanNormAtom) = ConvexVexity()

evaluate(x::EuclideanNormAtom) = norm(vec(evaluate(x.children[1])))

function new_conic_form!(context::Context{T}, A::EuclideanNormAtom) where {T}
    x = only(AbstractTrees.children(A))
    t = conic_form!(context, Variable())
    obj = conic_form!(context, x)
    f = operate(vcat, T, sign(A), t, obj)
    MOI_add_constraint(context.model, f, MOI.SecondOrderCone(length(x) + 1))
    return t
end

function LinearAlgebra.norm2(x::AbstractExpr)
    if sign(x) == ComplexSign()
        return EuclideanNormAtom([real(x); imag(x)])
    end
    return EuclideanNormAtom(x)
end
