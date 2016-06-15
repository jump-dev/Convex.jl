#############################################################################
# real_imag.jl
# Handles real and imaginary part of the variables, constants
# and expressions.
#############################################################################

import Base.real, Base.imag
export real, imag
export sign, monotonicity, curvature, evaluate, conic_form!


### Real
type RealAtom <: AbstractExpr
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
  return monotonicity(x.children[1])
end

function curvature(x::RealAtom)
  return ConstVexity()
end

function evaluate(x::RealAtom)
  return real(evaluate(x.children[1]))
end

function conic_form!(x::RealAtom, unique_conic_forms::UniqueConicForms)
  if !has_conic_form(unique_conic_forms, x)
    objective = conic_form!(x.children[1], unique_conic_forms)

    for var in keys(objective)
      x = real(objective[var][1])
      y = real(objective[var][2])
      objective[var] = (x,y)
    end
    cache_conic_form!(unique_conic_forms, x, objective)
  end
  return get_conic_form(unique_conic_forms, x)
end

real(x::AbstractExpr) = RealAtom(x)
real(x::Value) = RealAtom(Constant(x))



### Imaginary
type ImaginaryAtom <: AbstractExpr
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
  if sign(x.children[1]) == ComplexSign()
    return NoSign()
  else 
    return sign(Constant(0))
  end
end

function monotonicity(x::ImaginaryAtom)
  return monotonicity(x.children[1])
end

function curvature(x::ImaginaryAtom)
  return ConstVexity()
end

function evaluate(x::ImaginaryAtom)
  return imag(evaluate(x.children[1]))
end

function conic_form!(x::ImaginaryAtom, unique_conic_forms::UniqueConicForms)
  if !has_conic_form(unique_conic_forms, x)
    objective = conic_form!(x.children[1], unique_conic_forms)

    for var in keys(objective)
      x = imag(objective[var][1])
      y = imag(objective[var][2])
      objective[var] = (x,y)
    end
    cache_conic_form!(unique_conic_forms, x, objective)
  end
  return get_conic_form(unique_conic_forms, x)
end

imag(x::AbstractExpr) = ImaginaryAtom(x)
imag(x::Value) = ImaginaryAtom(Constant(x))
