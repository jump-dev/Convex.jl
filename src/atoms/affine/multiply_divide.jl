#############################################################################
# multiply_divide.jl
# Handles scalar multiplication, matrix multiplication, and scalar division
# of variables, constants and expressions.
# All expressions and atoms are subtpyes of AbstractExpr.
# Please read expressions.jl first.
#############################################################################

import Base.*, Base.(.*), Base./, Base.(./)
export *, .*
export sign, monotonicity, curvature, evaluate, conic_form!

### Scalar and matrix multiplication

type MultiplyAtom <: AbstractExpr
  head::Symbol
  id_hash::UInt64
  children::Tuple{AbstractExpr, AbstractExpr}
  size::Tuple{Int, Int}

  function MultiplyAtom(x::AbstractExpr, y::AbstractExpr)
    if x.size == (1, 1)
      sz = y.size
    elseif y.size == (1, 1)
      sz = x.size
    elseif x.size[2] ==  y.size[1]
      sz = (x.size[1], y.size[2])
    else
      error("Cannot multiply two expressions of sizes $(x.size) and $(y.size)")
    end
    children = (x, y)
    return new(:*, hash(children), children, sz)
  end
end

function sign(x::MultiplyAtom)
  return sign(x.children[1]) * sign(x.children[2])
end

function monotonicity(x::MultiplyAtom)
  return (sign(x.children[2]) * Nondecreasing(), sign(x.children[1]) * Nondecreasing())
end

# Multiplication has an indefinite hessian, so if neither children are constants,
# the curvature of the atom will violate DCP.
function curvature(x::MultiplyAtom)
  if vexity(x.children[1]) != ConstVexity() && vexity(x.children[2]) != ConstVexity()
    return NotDcp()
  else
    return ConstVexity()
  end
end

function evaluate(x::MultiplyAtom)
  return evaluate(x.children[1]) * evaluate(x.children[2])
end

function conic_form!(x::MultiplyAtom, unique_conic_forms::UniqueConicForms)
  if !has_conic_form(unique_conic_forms, x)
    # scalar multiplication
    if x.children[1].size == (1, 1) || x.children[2].size == (1, 1)
      if vexity(x.children[1]) == ConstVexity()
        const_child = x.children[1]
        expr_child = x.children[2]
      elseif vexity(x.children[2]) == ConstVexity()
        const_child = x.children[2]
        expr_child = x.children[1]
      else
        error("multiplication of two non-constant expressions is not DCP compliant")
      end
      objective = conic_form!(expr_child, unique_conic_forms)

      # make sure all 1x1 sized objects are interpreted as scalars, since
      # [1] * [1, 2, 3] is illegal in julia, but 1 * [1, 2, 3] is ok
      if const_child.size == (1, 1)
        const_multiplier = evaluate(const_child)[1]
      else
        const_multiplier = reshape(evaluate(const_child), get_vectorized_size(const_child), 1)
      end

      objective = const_multiplier * objective

    # left matrix multiplication
    elseif x.children[1].head == :constant
      objective = conic_form!(x.children[2], unique_conic_forms)
      objective = kron(speye(x.size[2]), x.children[1].value) * objective
    # right matrix multiplication
    else
      objective = conic_form!(x.children[1], unique_conic_forms)
      objective = kron(x.children[2].value', speye(x.size[1])) * objective
    end
    cache_conic_form!(unique_conic_forms, x, objective)
  end
  return get_conic_form(unique_conic_forms, x)
end

function *(x::AbstractExpr, y::AbstractExpr)
  if hash(x) == hash(y) && x.size == (1, 1)
    return square(x)
  end
  return MultiplyAtom(x, y)
end

*(x::Value, y::AbstractExpr) = MultiplyAtom(Constant(x), y)
*(x::AbstractExpr, y::Value) = MultiplyAtom(x, Constant(y))
/(x::AbstractExpr, y::Value) = MultiplyAtom(x, Constant(1./y))

### .*
# All constructors of this check (and so this function requires)
# that the first child be constant to have the expression be DCP
type DotMultiplyAtom <: AbstractExpr
  head::Symbol
  id_hash::UInt64
  children::Tuple{AbstractExpr, AbstractExpr}
  size::Tuple{Int, Int}

  function DotMultiplyAtom(x::AbstractExpr, y::AbstractExpr)
    # we don't explictly check, yet, that the sizes of x and y are compatible
    # if they're not, we'll get a DimensionMismatch error when we promote
    # the size of the coefficient matrix below
    # to check here, just run 
    # coeff = x.value .* ones(size(y))
    children = (x, y)
    return new(:.*, hash(children), children, y.size)
  end
end

function sign(x::DotMultiplyAtom)
  return sign(x.children[1]) * sign(x.children[2])
end

function monotonicity(x::DotMultiplyAtom)
  return (sign(x.children[2]) * Nondecreasing(), sign(x.children[1]) * Nondecreasing())
end

function curvature(x::DotMultiplyAtom)
  if vexity(x.children[1]) == ConstVexity()
    return ConstVexity()
  else
    return NotDcp()
  end
end

function evaluate(x::DotMultiplyAtom)
  return evaluate(x.children[1]) .* evaluate(x.children[2])
end

function conic_form!(x::DotMultiplyAtom, unique_conic_forms::UniqueConicForms)
  if !has_conic_form(unique_conic_forms, x)
    if vexity(x.children[1]) != ConstVexity()
      error("multiplication of two non-constant expressions is not DCP compliant")
    end    
    # promote the size of the coefficient matrix, so eg
    # 3 .* x 
    # works regardless of the size of x
    coeff = x.children[1].value .* ones(size(x.children[2]))
    const_multiplier = spdiagm(vec(coeff)) # used to be spdiagm((vec(coeff),), (0,)), not sure why
    objective = const_multiplier * conic_form!(x.children[2], unique_conic_forms)
    cache_conic_form!(unique_conic_forms, x, objective)
  end
  return get_conic_form(unique_conic_forms, x)
end

function .*(x::Constant, y::AbstractExpr)
  if x.size == (1, 1) || y.size == (1, 1)
    return x * y
  else
    return DotMultiplyAtom(x, y)
  end
end

.*(x::AbstractExpr, y::AbstractExpr) = vexity(x) == ConstVexity() ? DotMultiplyAtom(x, y) : DotMultiplyAtom(y, x)
.*(x::Value, y::AbstractExpr) = DotMultiplyAtom(Constant(x), y)
.*(x::AbstractExpr, y::Value) = DotMultiplyAtom(Constant(y), x)
./(x::AbstractExpr, y::Value) = DotMultiplyAtom(Constant(1./y), x)
