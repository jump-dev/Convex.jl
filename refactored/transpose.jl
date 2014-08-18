#############################################################################
# transpose.jl
# Returns the transpose of a matrix
# All expressions and atoms are subtpyes of AbstractExpr.
# Please read expressions.jl first.
#############################################################################

import Base.transpose, Base.ctranspose
export transpose, ctranspose, TransposeAtom
export sign, curvature, monotonicity, evaluate, dual_conic_form

type TransposeAtom <: AbstractExpr
  head::Symbol
  children_hash::Uint64
  children::(AbstractExpr,)
  size::(Int64, Int64)

  function TransposeAtom(x::AbstractExpr)
    children = (x,)
    return new(:transpose, hash(children), children, (x.size[2], x.size[1]))
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

# Since everything is vectorized, we simply need to multiply x by a permutation
# matrix such that coeff * vectorized(x) - vectorized(x') = 0
function dual_conic_form(x::TransposeAtom)
  objective, constraints = dual_conic_form(x.children[1])
  sz = get_vectorized_size(x.size)
  transpose_matrix = spzeros(sz, sz)
  num_rows = x.size[1]
  num_cols = x.size[2]

  for r = 1:num_rows
    for c = 1:num_cols
      i = (c - 1) * num_rows + r
      j = (r - 1) * num_cols + c
      transpose_matrix[i, j] = 1.0
    end
  end
  objective = transpose_matrix * objective
  return (objective, constraints)
end

transpose(x::AbstractExpr) = TransposeAtom(x)
ctranspose(x::AbstractExpr) = transpose(x)
