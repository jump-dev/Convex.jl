export geo_mean

function geo_mean(x::AbstractCvxExpr, y::AbstractCvxExpr)
  #TODO vexity and sign checks
  this = CvxExpr(:geo_mean, [x, y], :concave, :pos, x.size)
  x_size = get_vectorized_size(x)

  canon_constr_array = CanonicalConstr[]
  for i = 1:x_size
    coeffs1 = spzeros(3, x_size)
    coeffs1[1, i] = -1
    coeffs1[1, i] = 1
    coeffs2 = spzeros(3, x_size)
    coeffs2[1, i] = -1
    coeffs2[2, i] = -1
    coeffs3 = spzeros(3, x_size)
    coeffs3[3, i] = -2
    cone_coeffs = VecOrMatOrSparse[coeffs1, coeffs2, coeffs3]
    cone_vars = [x.uid, y.uid, this.uid]
    cone_constant = zeros(3, 1)
    push!(canon_constr_array, CanonicalConstr(cone_coeffs, cone_vars, cone_constant, false, true))
  end

  # x,y >= 0 linear constraint
  lin_coeffs1 = VecOrMatOrSparse[-speye(x_size)]
  lin_vars1 = [y.uid]
  lin_constant1 = zeros(x_size, 1)
  lin_coeffs2 = VecOrMatOrSparse[-speye(x_size)]
  lin_vars2 = [x.uid]
  lin_constant2 = zeros(x_size, 1)
  push!(canon_constr_array, CanonicalConstr(lin_coeffs1, lin_vars1, lin_constant1, false, false))
  push!(canon_constr_array, CanonicalConstr(lin_coeffs2, lin_vars2, lin_constant2, false, false))
  append!(canon_constr_array, x.canon_form())
  append!(canon_constr_array, y.canon_form())
  this.canon_form = ()->canon_constr_array
  this.evaluate = ()->Base.sqrt(x.evaluate() .* y.evaluate())

  return this
end
