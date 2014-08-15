export +, -
export sign, intrinsic_vexity, monotonicity, evaluate

### Unary Negation

type NegateAtom <: AbstractExpr
  head::Symbol
  id::Uint64
  children::(AbstractExpr,)
  size::(Int64, Int64)

  function NegateAtom(x::AbstractExpr)
    return new(:-, object_id(x), (x,), x.size)
  end
end

function sign(x::NegateAtom)
  return -sign(x.children[1])
end

function monotonicity(x::NegateAtom)
  return (Nonincreasing(),)
end

function intrinsic_vexity(x::NegateAtom)
  return ConstVexity()
end

function evaluate(x::NegateAtom)
  return -evaluate(x.children[1])
end

-(x::AbstractExpr) = NegateAtom(x)

function dual_conic_form(e::NegateAtom)
  objective, constraints = dual_conic_form(e.children[1])
  return (-objective, constraints)
end

### Binary Addition/Subtraction

type AdditionAtom <: AbstractExpr
  head::Symbol
  id::Uint64
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
    return new(:+, object_id(x) + object_id(y), (x, y), sz)
  end
end

function monotonicity(x::AdditionAtom)
  return (Nondecreasing(), Nondecreasing())
end

function intrinsic_vexity(x::AdditionAtom)
  return ConstVexity()
end

function evaluate(x::AdditionAtom)
  return evaluate(x.children[1]) + evaluate(x.children[2])
end


function dual_conic_form(x::AdditionAtom)
  child_cones = map(dual_conic_form, x.children)
  objective = Coefficients()
  constraints = ConeContr[]
  for (child_objective, child_constraints) in child_cones
    append!(constraints, child_constraints)
    objective += child_objective
  end
  return (objective, constraints)
end

+(x::AbstractExpr, y::AbstractExpr) = AdditionAtom(x, y)
+(x::Value, y::AbstractExpr) = AdditionAtom(Constant(x), y)
+(x::AbstractExpr, y::Value) = AdditionAtom(x, Constant(y))
-(x::AbstractExpr, y::AbstractExpr) = x + (-y)
-(x::Value, y::AbstractExpr) = Constant(x) + (-y)
-(x::AbstractExpr, y::Value) = x + (-Constant(y))
