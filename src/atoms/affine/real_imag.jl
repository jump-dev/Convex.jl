mutable struct RealAtom <: AbstractExpr
    children::Tuple{AbstractExpr}
    size::Tuple{Int,Int}

    RealAtom(x::AbstractExpr) = new((x,), x.size)
end

head(io::IO, ::RealAtom) = print(io, "real")

function Base.sign(x::RealAtom)
    if sign(x.children[1]) == ComplexSign()
        return NoSign()
    end
    return sign(x.children[1])
end

monotonicity(::RealAtom) = (Nondecreasing(),)

curvature(::RealAtom) = ConstVexity()

evaluate(x::RealAtom) = real.(evaluate(x.children[1]))

function new_conic_form!(context::Context{T}, x::RealAtom) where {T}
    obj = conic_form!(context, only(AbstractTrees.children(x)))
    return operate(real, T, sign(x), obj)
end

Base.real(x::AbstractExpr) = RealAtom(x)

Base.real(x::ComplexConstant) = x.real_constant

Base.real(x::Constant) = x

mutable struct ImaginaryAtom <: AbstractExpr
    children::Tuple{AbstractExpr}
    size::Tuple{Int,Int}

    ImaginaryAtom(x::AbstractExpr) = new((x,), x.size)
end

head(io::IO, ::ImaginaryAtom) = print(io, "imag")

Base.sign(::ImaginaryAtom) = NoSign()

monotonicity(::ImaginaryAtom) = (Nondecreasing(),)

curvature(::ImaginaryAtom) = ConstVexity()

evaluate(x::ImaginaryAtom) = imag.(evaluate(x.children[1]))

function new_conic_form!(context::Context{T}, x::ImaginaryAtom) where {T}
    obj = conic_form!(context, only(AbstractTrees.children(x)))
    return operate(imag, T, sign(x), obj)
end

Base.imag(x::AbstractExpr) = ImaginaryAtom(x)

Base.imag(x::ComplexConstant) = x.imag_constant

Base.imag(x::Constant) = Constant(zero(x.value))
