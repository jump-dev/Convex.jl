import Base.transpose, Base.ctranspose
export transpose, ctranspose

function transpose(x::AbstractCvxExpr)
  sz = get_vectorized_size(x.size)
  coeffs = spzeros(sz, sz)
  num_rows = x.size[1]
  num_cols = x.size[2]

  for r = 1:num_rows
    for c = 1:num_cols
      i = (c - 1) * num_rows + r
      j = (r - 1) * num_cols + c
      coeffs[i, j] = 1.0
    end
  end

  this = CvxExpr(:transpose, [x], x.vexity, x.sign, (x.size[2], x.size[1]))

  coeffs = VecOrMatOrSparse[speye(sz), -coeffs]
  vars = [x.uid, this.uid]
  canon_constr = CanonicalConstr(coeffs, vars, spzeros(sz, 1), true, false)

  canon_constr_array = x.canon_form()
  push!(canon_constr_array, canon_constr)

  this.canon_form = ()->canon_constr_array
  this.evaluate = ()->x.evaluate()'
  return this
end

ctranspose(x::AbstractCvxExpr) = transpose(x)

function transpose(x::Constant)
  return Constant(x.value')
end

ctranspose(x::Constant) = transpose(x)
