import Base.diag
export diag

# Finds the "k"-th diagonal of x as a column vector.
# If k == 0, it returns the main diagonal and so on.
# Let x be of size m x n and d be the diagonal.
# Since x is vectorized, the way canonicalization works is:
#
# 1. We calculate the size of the diagonal (sz_diag) and the first index
# of vectorized x that will be part of d
# 2. We create the coefficient matrix for vectorized x, called coeff of size
# sz_diag x mn
# 3. We populate coeff with 1s at the correct indices
# The canonical form will then be:
# coeff * x - d = 0
function diag(x::AbstractCvxExpr, k::Int64=0)
  (num_rows, num_cols) = x.size

  if k >= num_cols || k <= -num_rows
    error("Bounds error in calling diag")
  end

  if k >= 0
    start_index = k * num_rows + 1
    sz_diag = Base.min(num_rows, num_cols - k)
  else
    start_index = -k + 1
    sz_diag = Base.min(num_rows + k, num_cols)
  end

  coeff = spzeros(sz_diag, get_vectorized_size(x))

  for i in 1:sz_diag
    coeff[i, start_index] = 1
    start_index += num_rows + 1
  end

  this = CvxExpr(:diag, [x], x.vexity, x.sign, (sz_diag, 1))
  coeffs = VecOrMatOrSparse[coeff, -speye(sz_diag)]
  vars = [x.uid, this.uid]
  constant = zeros(sz_diag, 1)

  canon_constr_array = [CanonicalConstr(coeffs, vars, constant, true, false)]
  append!(canon_constr_array, x.canon_form())
  this.canon_form = ()->canon_constr_array

  this.evaluate = ()->Base.diag(x.evaluate(), k)
  return this
end
