mutable struct ConjugateAtom <: AbstractExpr
    children::Tuple{AbstractExpr}
    size::Tuple{Int,Int}

    ConjugateAtom(x::AbstractExpr) = new((x,), x.size)
end

head(io::IO, ::ConjugateAtom) = print(io, "conj")

Base.sign(x::ConjugateAtom) = sign(x.children[1])

monotonicity(::ConjugateAtom) = (Nondecreasing(),)

curvature(::ConjugateAtom) = ConstVexity()

evaluate(x::ConjugateAtom) = conj(evaluate(x.children[1]))

function new_conic_form!(context::Context{T}, x::ConjugateAtom) where {T}
    objective = conic_form!(context, only(AbstractTrees.children(x)))
    return operate(conj, T, sign(x), objective)
end

function Base.conj(x::AbstractExpr)
    if sign(x) == ComplexSign()
        return ConjugateAtom(x)
    end
    return x
end

Base.conj(x::Constant) = x

Base.conj(x::ComplexConstant) = ComplexConstant(real(x), -imag(x))
