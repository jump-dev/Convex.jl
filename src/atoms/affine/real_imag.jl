#############################################################################
# real_imag.jl
# Handles real and imaginary part of the variables, constants
# and expressions.
#############################################################################

import Base.real, Base.imag
export real, imag
export sign, curvature, monotonicity, evaluate


### Real
type RealAtom <: AbstractExpr
  head::Symbol
  id_hash::UInt64
  children::Array{AbstractExpr, 1}
  size::Tuple{Int, Int}

  function RealAtom(x::AbstractExpr)
    children = (x,)
    return new(:real, hash(children), children, x.size)
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
    objective = ConicObj()
    for child in x.children
      child_objective = conic_form!(child, unique_conic_forms)
      if x.size != child.size
        child_objective = promote_size(child_objective, get_vectorized_size(x))
      end
      objective += child_objective
    end
    cache_conic_form!(unique_conic_forms, x, objective)
  end
  return get_conic_form(unique_conic_forms, x)
end

real(x::AbstractExpr) = RealAtom(x)
real(x::Value) = RealAtom(Constant(x))
