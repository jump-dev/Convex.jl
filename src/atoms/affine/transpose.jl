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
  canon_constr_array = Any[{
    :coeffs => Any[speye(sz), -coeffs],
    :vars => [x.uid, this.uid],
    :constant => spzeros(sz, 1),
    :is_eq => true
  }]

  append!(canon_constr_array, x.canon_form())
  this.canon_form = ()->canon_constr_array
  return this
end

ctranspose(x::AbstractCvxExpr) = transpose(x)

function transpose(x::Constant)
  return Constant(x.value')
end

ctranspose(x::Constant) = transpose(x)
