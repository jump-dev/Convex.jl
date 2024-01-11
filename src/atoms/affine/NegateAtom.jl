mutable struct NegateAtom <: AbstractExpr
    children::Tuple{AbstractExpr}
    size::Tuple{Int,Int}

    NegateAtom(x::AbstractExpr) = new((x,), x.size)
end

Base.sign(x::NegateAtom) = -sign(x.children[1])

monotonicity(::NegateAtom) = (Nonincreasing(),)

curvature(::NegateAtom) = ConstVexity()

evaluate(x::NegateAtom) = -evaluate(x.children[1])

function new_conic_form!(context::Context{T}, A::NegateAtom) where {T}
    subobj = conic_form!(context, only(AbstractTrees.children(A)))
    if subobj isa Value
        return -subobj
    end
    return operate(-, T, sign(A), subobj)
end

Base.:-(x::AbstractExpr) = NegateAtom(x)

Base.:-(x::Union{Constant,ComplexConstant}) = constant(-evaluate(x))
