#############################################################################
# add_subtract.jl
# Handles unary negation, addition and subtraction of variables, constants
# and expressions.
# All expressions and atoms are subtpyes of AbstractExpr.
# Please read expressions.jl first.
#############################################################################

export +, -
export sign, curvature, monotonicity, evaluate

### Unary Negation

type NegateAtom <: AbstractExpr
  head::Symbol
  children_hash::Uint64
  children::(AbstractExpr,)
  size::(Int64, Int64)

  function NegateAtom(x::AbstractExpr)
    children = (x,)
    return new(:-, hash(children), children, x.size)
  end
end

function sign(x::NegateAtom)
  return -sign(x.children[1])
end

function monotonicity(x::NegateAtom)
  return (Nonincreasing(),)
end

function curvature(x::NegateAtom)
  return ConstVexity()
end

function evaluate(x::NegateAtom)
  return -evaluate(x.children[1])
end

-(x::AbstractExpr) = NegateAtom(x)

function conic_form(x::NegateAtom, unique_constr)
  if !((x.head, x.children_hash) in unique_constr)
    objective, constraints = conic_form(x.children[1], unique_constr)
    objective = -objective
    unique_constr[(x.head, x.children_hash)] = (objective, constraints)
  end
  return unique_constr[(x.head, x.children_hash)]
end


### Binary Addition/Subtraction
type AdditionAtom <: AbstractExpr
  head::Symbol
  children_hash::Uint64
  children::(AbstractExpr, AbstractExpr)
  size::(Int64, Int64)

  function AdditionAtom(x::AbstractExpr, y::AbstractExpr)
    if x.size == y.size || y.size == (1, 1)
      sz = x.size
    elseif x.size == (1, 1)
      sz = y.size
    else
      error("Cannot add expressions of sizes $(x.size) and $(y.size)")
    end
    children = (x, y)
    return new(:+, hash(children), children, sz)
  end
end

function sign(x::AdditionAtom)
  return sign(x.children[1]) + sign(x.children[2])
end

function monotonicity(x::AdditionAtom)
  return (Nondecreasing(), Nondecreasing())
end

function curvature(x::AdditionAtom)
  return ConstVexity()
end

function evaluate(x::AdditionAtom)
  return evaluate(x.children[1]) + evaluate(x.children[2])
end

function conic_form(x::AdditionAtom, unique_constr)
  if !((x.head, x.children_hash) in unique_constr)
    objective, constraints = conic_form(x.children[1], unique_constr)
    objective2, constraints2 = conic_form(x.children[2], unique_constr)
    append!(constraints, constraints2)
    objective += objective2
    unique_constr[(x.head, x.children_hash)] = (objective, constraints)
  end
  return unique_constr[(x.head, x.children_hash)]
end

+(x::AbstractExpr, y::AbstractExpr) = AdditionAtom(x, y)
+(x::Value, y::AbstractExpr) = AdditionAtom(Constant(x), y)
+(x::AbstractExpr, y::Value) = AdditionAtom(x, Constant(y))
-(x::AbstractExpr, y::AbstractExpr) = x + (-y)
-(x::Value, y::AbstractExpr) = Constant(x) + (-y)
-(x::AbstractExpr, y::Value) = x + Constant(-y)
