export *
export sign, monotonicity, intrinsic_vexity

### Multiplication

type MutiplyAtom <: AbstractExpr
  head::Symbol
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
    return new(:*, (x, y), sz)
  end
end

function sign(x::MultiplyAtom)
  return sign(x.children[1]) * sign(x.children[2])
end

function monotonicity(x::MultiplyAtom)
  return (sign(x.children[2]) * Nondecreasing(), sign(x.children[1]) * Nondecreasing())
end

function intrinsic_vexity(x::MultiplyAtom)
  if vexity(x.children[1]) != Constant() && vexity(x.children[2]) != Constant()
    return NoVexity()
  else
    return Affine()
  end
end
