struct ConjugateAtom <: AbstractExpr
    head::Symbol
    id_hash::UInt64
    children::Tuple{AbstractExpr}
    size::Tuple{Int,Int}

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

function conic_form!(context::Context{T}, x::ConjugateAtom) where {T}
    objective = conic_form!(context, only(children(x)))
    return operate(conj, T, sign(x), objective)
end

function Base.conj(x::AbstractExpr)
    if sign(x) == ComplexSign()
        return ConjugateAtom(x)
    else
        return x
    end
end
Base.conj(x::Constant) = x
Base.conj(x::ComplexConstant) = ComplexConstant(real(x), -imag(x))
