#############################################################################
# real_imag.jl
# Handles real and imaginary part of the variables, constants
# and expressions.
#############################################################################

### Real
mutable struct RealAtom <: AbstractExpr
    children::Tuple{AbstractExpr}
    size::Tuple{Int,Int}

    function RealAtom(x::AbstractExpr)
        children = (x,)
        return new(children, x.size)
    end
end

head(io::IO, ::RealAtom) = print(io, "real")

function Base.sign(x::RealAtom)
    if sign(x.children[1]) == ComplexSign()
        return NoSign()
    else
        return sign(x.children[1])
    end
end

function monotonicity(x::RealAtom)
    return (Nondecreasing(),)
end

function curvature(x::RealAtom)
    return ConstVexity()
end

function evaluate(x::RealAtom)
    return real.(evaluate(x.children[1]))
end

function new_conic_form!(context::Context{T}, x::RealAtom) where {T}
    obj = conic_form!(context, only(AbstractTrees.children(x)))
    return operate(real, T, sign(x), obj)
end

Base.real(x::AbstractExpr) = RealAtom(x)
Base.real(x::ComplexConstant) = x.real_constant
Base.real(x::Constant) = x

### Imaginary
mutable struct ImaginaryAtom <: AbstractExpr
    children::Tuple{AbstractExpr}
    size::Tuple{Int,Int}

    function ImaginaryAtom(x::AbstractExpr)
        children = (x,)
        return new(children, x.size)
    end
end

head(io::IO, ::ImaginaryAtom) = print(io, "imag")

function Base.sign(x::ImaginaryAtom)
    sign(x.children[1]) == ComplexSign()
    return NoSign()
end

function monotonicity(x::ImaginaryAtom)
    return (Nondecreasing(),)
end

function curvature(x::ImaginaryAtom)
    return ConstVexity()
end

function evaluate(x::ImaginaryAtom)
    return imag.(evaluate(x.children[1]))
end

function new_conic_form!(context::Context{T}, x::ImaginaryAtom) where {T}
    obj = conic_form!(context, only(AbstractTrees.children(x)))
    return operate(imag, T, sign(x), obj)
end

Base.imag(x::AbstractExpr) = ImaginaryAtom(x)
Base.imag(x::ComplexConstant) = x.imag_constant
Base.imag(x::Constant) = Constant(zero(x.value))
