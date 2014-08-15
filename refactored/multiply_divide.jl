export *
export sign, monotonicity, intrinsic_vexity, evaluate

### Multiplication

type MultiplyAtom <: AbstractExpr
  head::Symbol
  id::Uint64
  children::(AbstractExpr, AbstractExpr)
  size::(Int64, Int64)

  function MutiplyAtom(x::AbstractExpr, y::AbstractExpr)
    if x.size == (1, 1)
      sz = y.size
    elseif y.size == (1, 1)
      sz = x.size
    elseif x.size[2] ==  y.size[1]
      sz = (x.size[1], y.size[2])
    else
      error("Cannot multiply two expressions of sizes $(x.size) and $(y.size)")
    end
    return new(:*, object_id(x) + object_id(y), (x, y), sz)
  end
end

function sign(x::MultiplyAtom)
  return sign(x.children[1]) * sign(x.children[2])
end

function monotonicity(x::MultiplyAtom)
  return (sign(x.children[2]) * Nondecreasing(), sign(x.children[1]) * Nondecreasing())
end

function intrinsic_vexity(x::MultiplyAtom)
  if vexity(x.children[1]) != ConstVexity() && vexity(x.children[2]) != ConstVexity()
    return NoVexity()
  else
    return ConstVexity()
  end
end

function evaluate(x::MultiplyAtom)
  return evaluate(x.children[1]) * evaluate(x.children[2])
end

function dual_conic_form(x::MultiplyAtom)
end
