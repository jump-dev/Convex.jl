export +, -
export sign, vexity, intrinsic_vexity, monotonicity, evaluate

### Unary Negation

type NegateAtom <: AffineFunc
  head::Symbol
  children::(AbstractExpr,)
  size::(Int64, Int64)

  function NegateAtom(x::AbstractExpr)
    return new(:-, (x,), x.size)
  end
end

function sign(x::NegateAtom)
  return -sign(x.children[1])
end

function vexity(x::NegateAtom)
  return -vexity(x.children[1])
end

function monotonicity(x::NegateAtom)
  return (Nonincreasing(),)
end

function evaluate(x::NegateAtom)
  return -evaluate(x.children[1])
end

-(x::AbstractExpr) = NegateAtom(x)


### Binary Addition/Subtraction

type AdditionAtom <: AffineFunc
  head::Symbol
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
    return new(:+, (x, y), sz)
end

function monotonicity(x::AdditionAtom)
  return (Nondecreasing(), Nondecreasing())
end

function intrinsic_vexity(x::AdditionAtom)
  return Affine()
end

function evaluate(x::AdditionAtom)
  return evaluate(x.children[1]) + evaluate(x.children[2])
end

+(x::AbstractExpr, y::AbstractExpr) = AdditionAtom(x, y)
+(x::Value, y::AbstractExpr) = AdditionAtom(Constant(x), y)
+(x::AbstractExpr, y::Value) = AdditionAtom(x, Constant(y))
-(x::AbstractExpr, y::AbstractExpr) = x + (-y)
-(x::Value, y::AbstractExpr) = Constant(x) + (-y)
-(x::AbstractExpr, y::Value) = x + (-Constant(y))
