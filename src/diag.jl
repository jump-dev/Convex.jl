#############################################################################
# diag.jl
# Returns the kth diagonal of a matrix expression
# All expressions and atoms are subtpyes of AbstractExpr.
# Please read expressions.jl first.
#############################################################################

import Base.diag
export diag

### Diagonal
### Represents the kth diagonal of an mxn matrix as a (min(m, n) - k) x 1 vector
type DiagAtom <: AbstractExpr
  head::Symbol
  children_hash::Uint64
  children::(AbstractExpr, Int64)
  size::(Int64, Int64)

  function DiagAtom(x::AbstractExpr, k::Int64=0)
    (num_rows, num_cols) = x.size

    if k >= num_cols || k <= -num_rows
      error("Bounds error in calling diag")
    end

    children = (x, k)
    return new(:sum, hash(children), children, (minimum(x.size) - k, 1))
  end
end

function sign(x::DiagAtom)
  return sign(x.children[1])
end

# The monotonicity
function monotonicity(x::DiagAtom)
  return (Nondecreasing(),)
end

# If we have h(x) = f o g(x), the chain rule says h''(x) = g'(x)^T f''(g(x))g'(x) + f'(g(x))g''(x);
# this represents the first term
function curvature(x::DiagAtom)
  return ConstVexity()
end

function evaluate(x::DiagAtom)
  return diag(evaluate(x.children[1]), x.children[2])
end

diag(x::AbstractExpr, k::Int64=0) = DiagAtom(x, k)

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
function conic_form(x::DiagAtom, unique_constr)
  if !((x.head, x.children_hash) in unique_constr)
    (num_rows, num_cols) = x.children[1].size
    k = x.children[2]

    if k >= 0
      start_index = k * num_rows + 1
      sz_diag = Base.min(num_rows, num_cols - k)
    else
      start_index = -k + 1
      sz_diag = Base.min(num_rows + k, num_cols)
    end

    select_diag = spzeros(sz_diag, get_vectorized_size(x.children[1]))
    for i in 1:sz_diag
      select_diag[i, start_index] = 1
      start_index += num_rows + 1
    end

    objective, constraints = conic_form(x.children[1], unique_constr)
    new_obj = select_diag * objective
    unique_constr[(x.head, x.children_hash)] = (new_obj, constraints)
  end
  return unique_constr[(x.head, x.children_hash)]
end
