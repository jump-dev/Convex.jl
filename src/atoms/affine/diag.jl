import Base.diag
export diag

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
