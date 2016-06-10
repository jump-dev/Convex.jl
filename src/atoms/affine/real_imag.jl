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
    # see if we're forming a sum of more than two terms and condense them
    children = AbstractExpr[]
    if isa(x, AdditionAtom)
      append!(children, x.children)
    else
      push!(children, x)
    end
    if isa(y, AdditionAtom)
      append!(children, y.children)
    else
      push!(children, y)
    end
    return new(:+, hash(children), children, sz)
  end
end

function sign(x::RealAtom)
  return sum(Sign[sign(child) for child in x.children])
  # Creating an array of type Sign and adding all the sign of xhildren of x so if anyone is complex the resultant sign would be complex.
end

function monotonicity(x::RealAtom)
  return Monotonicity[Nondecreasing() for child in x.children]
end

function curvature(x::RealAtom)
  return ConstVexity()
end

function evaluate(x::RealAtom)
  return sum([evaluate(child) for child in x.children])
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
