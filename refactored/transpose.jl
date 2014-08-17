import Base.transpose, Base.ctranspose
export transpose, ctranspose, TransposeAtom
export sign, curvature, monotonicity, evaluate

# Since everything is vectorized, the canonical form of x' is simply
# multiplying x by a permutation matrix such that coeff * vectorized(x) - vectorized(x') = 0
type TransposeAtom <: AbstractExpr
  head::Symbol
  id::Uint64
  children::(AbstractExpr,)
  size::(Int64, Int64)

  function TransposeAtom(x::AbstractExpr)
    return new(:transpose, object_id(x), (x,), (x.size[2], x.size[1]))
  end
end

function sign(x::TransposeAtom)
  return sign(x.children[1])
end

function monotonicity(x::TransposeAtom)
  return (Nondecreasing(),)
end

function curvature(x::TransposeAtom)
  return ConstVexity()
end

function evaluate(x::TransposeAtom)
  return evaluate(x.children[1])'
end


transpose(x::AbstractExpr) = TransposeAtom(x)
ctranspose(x::AbstractExpr) = transpose(x)
