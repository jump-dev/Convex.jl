#############################################################################
# diagm.jl
# Converts a vector of size n into an n x n diagonal
# All expressions and atoms are subtpyes of AbstractExpr.
# Please read expressions.jl first.
#############################################################################

import Base.diagm
export diagm

type DiagMatrixAtom <: AbstractExpr
  head::Symbol
  id_hash::UInt64
  children::@compat Tuple{AbstractExpr}
  size::@compat Tuple{Int, Int}

  function DiagMatrixAtom(x::AbstractExpr)
    (num_rows, num_cols) = x.size

    if num_rows == 1
      sz = num_cols
    elseif num_cols == 1
      sz = num_rows
    else
      error("Only vectors are allowed for diagm. Did you mean to use diag?")
    end

    children = (x, )
    return new(:diagm, hash(children), children, (sz, sz))
  end
end

function sign(x::DiagMatrixAtom)
  return sign(x.children[1])
end

# The monotonicity
function monotonicity(x::DiagMatrixAtom)
  return (Nondecreasing(),)
end

# If we have h(x) = f o g(x), the chain rule says h''(x) = g'(x)^T f''(g(x))g'(x) + f'(g(x))g''(x);
# this represents the first term
function curvature(x::DiagMatrixAtom)
  return ConstVexity()
end

function evaluate(x::DiagMatrixAtom)
  return diagm(vec(evaluate(x.children[1])))
end

diagm(x::AbstractExpr) = DiagMatrixAtom(x)

function conic_form!(x::DiagMatrixAtom, unique_conic_forms::UniqueConicForms)
  if !has_conic_form(unique_conic_forms, x)
    sz = x.size[1]

    I = 1:sz+1:sz*sz
    J = 1:sz
    coeff = sparse(I, J, 1.0, sz * sz, sz)
    objective = conic_form!(x.children[1], unique_conic_forms)
    new_obj = coeff * objective
    cache_conic_form!(unique_conic_forms, x, new_obj)
  end
  return get_conic_form(unique_conic_forms, x)
end
