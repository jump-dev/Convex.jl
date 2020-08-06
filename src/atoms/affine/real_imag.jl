#############################################################################
# real_imag.jl
# Handles real and imaginary part of the variables, constants
# and expressions.
#############################################################################

import Base.real, Base.imag

### Real
struct RealAtom <: AbstractExpr
    head::Symbol
    id_hash::UInt64
    children::Tuple{AbstractExpr}
    size::Tuple{Int, Int}

    function RealAtom(x::AbstractExpr)
        children = (x,)
        return new(:real, hash(children), children, x.size)
    end
end

function sign(x::RealAtom)
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

function template(x::RealAtom, context::Context{T}) where T
    obj = template(only(children(x)), context)
    return operate(real, T, obj)
end

real(x::AbstractExpr) = RealAtom(x)
real(x::Value) = Constant(real(x))
real(x::ComplexConstant) = x.real_constant
real(x::Constant) = x
real(x::ComplexVariable) = x.real_var



### Imaginary
struct ImaginaryAtom <: AbstractExpr
    head::Symbol
    id_hash::UInt64
    children::Tuple{AbstractExpr}
    size::Tuple{Int, Int}

    function ImaginaryAtom(x::AbstractExpr)
        children = (x,)
        return new(:imag, hash(children), children, x.size)
    end
end

function sign(x::ImaginaryAtom)
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

function template(x::ImaginaryAtom, context::Context{T}) where T
    obj = template(only(children(x)), context)
    return operate(imag, T, obj)
end

imag(x::AbstractExpr) = ImaginaryAtom(x)
imag(x::Value) = Constant(imag(x)) # is this `Constant` needed?
imag(x::ComplexVariable) = x.imag_var
imag(x::ComplexConstant) = x.imag_constant
imag(x::Constant) = Constant(zero(x.value))
