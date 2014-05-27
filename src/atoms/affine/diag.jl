export diag

# TODO: implement diag for non-main diagonal
function diag(x::AbstractCvxExpr)
  (num_rows, num_cols) = x.size

  num_rows_coeff = Base.min(num_rows, num_cols)
  num_cols_coeff = get_vectorized_size(x)

  coeff = spzeros(num_rows_coeff, num_cols_coeff)

  for i in num_rows_coeff
    idx = num_rows * (i - 1) + i
    coeff[i, idx] = 1
  end

  this = CvxExpr(:diag, [x], x.vexity, x.sign, (num_rows_coeff, 1))
  coeffs = VecOrMatOrSparse[coeff, -speye(num_rows_coeff)]
  vars = [x.uid, this.uid]
  constant = zeros(num_rows_coeff)

  canon_constr_array = [CanonicalConstr(coeffs, vars, constant, true, false)]
  append!(canon_constr_array, x.canon_form())
  this.canon_form = ()->canon_constr_array

  this.evaluate = ()->Base.diag(x.evaluate(), 0)
  return this
end
