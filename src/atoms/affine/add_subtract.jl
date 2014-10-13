#############################################################################
# add_subtract.jl
# Handles unary negation, addition and subtraction of variables, constants
# and expressions.
# All expressions and atoms are subtpyes of AbstractExpr.
# Please read expressions.jl first.
#############################################################################

export +, -, .+, .-
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
  if !((x.head, x.children_hash) in keys(unique_constr))
    objective, constraints = conic_form(x.children[1], unique_constr)
    objective = -objective
    unique_constr[(x.head, x.children_hash)] = (objective, constraints)
  end
  return safe_copy(unique_constr[(x.head, x.children_hash)])
end


### Addition
type AdditionAtom <: AbstractExpr
  head::Symbol
  children_hash::Uint64
  children::Array{AbstractExpr, 1}
  size::(Int64, Int64)

  function AdditionAtom(x::AbstractExpr, y::AbstractExpr)
    # find the size of the expression = max of size of x and size of y
    if x.size == y.size || y.size == (1, 1)
      sz = x.size
    elseif x.size == (1, 1)
      sz = y.size
    else
      error("Cannot add expressions of sizes $(x.size) and $(y.size)")
    end
    # see if we're forming a sum of more than two terms and condense them
    children = AbstractExpr[]
    if isa(x, AdditionAtom)
      append!(children, x.children)
      push!(children, y)
    elseif isa(y, AdditionAtom)
      append!(children, y.children)
      push!(children, x)
    else
      children = [x, y]
    end
    return new(:+, hash(children), children, sz)
  end
end

function sign(x::AdditionAtom)
  return sum(Sign[sign(child) for child in x.children])
end

function monotonicity(x::AdditionAtom)
  return Monotonicity[Nondecreasing() for child in x.children]
end

function curvature(x::AdditionAtom)
  return ConstVexity()
end

function evaluate(x::AdditionAtom)
  return sum([evaluate(child) for child in x.children])
end

function conic_form(x::AdditionAtom, unique_constr)
  if !((x.head, x.children_hash) in keys(unique_constr))
    objective, constraints = ConicObj(), ConicConstr[]
    for child in x.children
      child_objective, child_constraints = conic_form(child, unique_constr)
      if x.size != child.size
        child_objective = promote_size(child_objective, get_vectorized_size(x))
      end
      append!(constraints, child_constraints)
      objective += child_objective
    end
    unique_constr[(x.head, x.children_hash)] = (objective, constraints)
  end
  return safe_copy(unique_constr[(x.head, x.children_hash)])
end

+(x::AbstractExpr, y::AbstractExpr) = AdditionAtom(x, y)
+(x::Value, y::AbstractExpr) = AdditionAtom(Constant(x), y)
+(x::AbstractExpr, y::Value) = AdditionAtom(x, Constant(y))
-(x::AbstractExpr, y::AbstractExpr) = x + (-y)
-(x::Value, y::AbstractExpr) = Constant(x) + (-y)
-(x::AbstractExpr, y::Value) = x + Constant(-y)

.+(x::AbstractExpr, y::AbstractExpr) = AdditionAtom(x, y)
.+(x::Value, y::AbstractExpr) = AdditionAtom(Constant(x), y)
.+(x::AbstractExpr, y::Value) = AdditionAtom(x, Constant(y))
.-(x::AbstractExpr, y::AbstractExpr) = x + (-y)
.-(x::Value, y::AbstractExpr) = Constant(x) + (-y)
.-(x::AbstractExpr, y::Value) = x + Constant(-y)
