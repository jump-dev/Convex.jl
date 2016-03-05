# How to create your own atom
# This file is a good example of how a user can create his or her own atom
# This antidiag atom is then included in n_queens.jl to solve the n queens problem

#############################################################################
# antidiag.jl
# Returns the kth anti-diagonal (counterdiagonal, secondary diagonal, or minor diagonal) of a matrix expression
# All expressions and atoms are subtpyes of AbstractExpr.
# Please read expressions.jl first.
#############################################################################

import Convex.sign, Convex.monotonicity, Convex.curvature, Convex.evaluate, Convex.conic_form!
export antidiag, sign, monotonicity, curvature, evaluate, conic_form!

### Diagonal
### Represents the kth diagonal of an mxn matrix as a (min(m, n) - k) x 1 vector
type AntidiagAtom <: AbstractExpr
  head::Symbol
  id_hash::UInt64
  children::Tuple{AbstractExpr}
  size::Tuple{Int, Int}
  k::Int

  function AntidiagAtom(x::AbstractExpr, k::Int=0)
    (num_rows, num_cols) = x.size

    if k >= num_cols || k <= -num_rows
      error("Bounds error in calling diag")
    end

    children = (x, )
    return new(:antidiag, hash((children, k)), children, (minimum(x.size) - k, 1), k)
  end
end

function sign(x::AntidiagAtom)
  return sign(x.children[1])
end

# The monotonicity
function monotonicity(x::AntidiagAtom)
  return (Nondecreasing(),)
end

# If we have h(x) = f o g(x), the chain rule says h''(x) = g'(x)^T f''(g(x))g'(x) + f'(g(x))g''(x);
# this represents the first term
function curvature(x::AntidiagAtom)
  return ConstVexity()
end

function evaluate(x::AntidiagAtom)
  return diag(flipdim(evaluate(x.children[1]),1), x.k)
end

antidiag(x::AbstractExpr, k::Int=0) = AntidiagAtom(x, k)

# Finds the "k"-th diagonal of x as a column vector
# If k == 0, it returns the main diagonal and so on
# Let x be of size m x n and d be the diagonal
# Since x is vectorized, the way canonicalization works is:
#
# 1. We calculate the size of the diagonal (sz_diag) and the first index
# of vectorized x that will be part of d
# 2. We create the coefficient matrix for vectorized x, called coeff of size
# sz_diag x mn
# 3. We populate coeff with 1s at the correct indices
# The canonical form will then be:
# coeff * x - d = 0
function conic_form!(x::AntidiagAtom, unique_conic_forms::UniqueConicForms)
  if !has_conic_form(unique_conic_forms, x)
    (num_rows, num_cols) = x.children[1].size
    k = x.k

    if k >= 0
      start_index = k * num_rows + num_rows
      sz_diag = Base.min(num_rows, num_cols - k)
    else
      start_index = num_rows + k
      sz_diag = Base.min(num_rows + k, num_cols)
    end

    select_diag = spzeros(sz_diag, get_vectorized_size(x.children[1]))
    for i in 1:sz_diag
      select_diag[i, start_index] = 1
      start_index += num_rows - 1
    end

    objective = conic_form!(x.children[1], unique_conic_forms)
    new_obj = select_diag * objective
    cache_conic_form!(unique_conic_forms, x, new_obj)
  end
  return get_conic_form(unique_conic_forms, x)
end
