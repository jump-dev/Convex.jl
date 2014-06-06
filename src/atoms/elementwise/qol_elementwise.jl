export qol_elementwise

function qol_elementwise(x::Constant, y::Constant)
  return Constant((x.value.^2)./y.value, :pos)
end

function qol_elementwise(x::Constant, y::AbstractCvxExpr)
  # TODO: vexity and sign checks
  this = CvxExpr(:qol_elementwise, [x, y], :convex, :pos, x.size)
  x_size = get_vectorized_size(x)

  canon_constr_array = CanonicalConstr[]
  for i = 1:x_size
    coeffs1 = spzeros(3, x_size)
    coeffs1[1, i] = -1
    coeffs1[2, i] = -1
    coeffs2 = spzeros(3, x_size)
    coeffs2[1, i] = -1
    coeffs2[2, i] = 1
    cone_coeffs = VecOrMatOrSparse[coeffs1, coeffs2]
    cone_vars = [y.uid, this.uid]
    cone_constant = [0; 0; 2*x.value[i]]
    push!(canon_constr_array, CanonicalConstr(cone_coeffs, cone_vars, cone_constant, false, true))
  end

  # y >= 0 linear constraint
  lin_coeffs = VecOrMatOrSparse[-speye(x_size)]
  lin_vars = [y.uid]
  lin_constant = zeros(x_size, 1)
  push!(canon_constr_array, CanonicalConstr(lin_coeffs, lin_vars, lin_constant, false, false))
  append!(canon_constr_array, y.canon_form())
  this.canon_form = ()->canon_constr_array
  this.evaluate = ()->x.evaluate() .^ 2 ./ y.evaluate()

  return this
end

function qol_elementwise(x::AbstractCvxExpr, y::Constant)
  # TODO: vexity and sign checks
  this = CvxExpr(:qol_elementwise, [x, y], :convex, :pos, x.size)
  x_size = get_vectorized_size(x)

  canon_constr_array = CanonicalConstr[]
  for i = 1:x_size
    coeffs1 = spzeros(3, x_size)
    coeffs1[3, i] = -2
    coeffs2 = spzeros(3, x_size)
    coeffs2[1, i] = -1
    coeffs2[2, i] = 1
    cone_coeffs = VecOrMatOrSparse[coeffs1, coeffs2]
    cone_vars = [x.uid, this.uid]
    cone_constant = [y.value[i], y.value[i], 0]
    push!(canon_constr_array, CanonicalConstr(cone_coeffs, cone_vars, cone_constant, false, true))
  end

  append!(canon_constr_array, x.canon_form())
  this.canon_form = ()->canon_constr_array
  this.evaluate = ()->x.evaluate() .^2 ./ y.evaluate()

  return this
end

function qol_elementwise(x::AbstractCvxExpr, y::AbstractCvxExpr)
  # TODO: vexity and sign checks
  this = CvxExpr(:qol_elementwise, [x, y], :convex, :pos, x.size)
  x_size = get_vectorized_size(x)

  canon_constr_array = CanonicalConstr[]
  for i = 1:x_size
    coeffs1 = spzeros(3, x_size)
    coeffs1[3, i] = -2
    coeffs2 = spzeros(3, x_size)
    coeffs2[1, i] = -1
    coeffs2[2, i] = -1
    coeffs3 = spzeros(3, x_size)
    coeffs3[1, i] = -1
    coeffs3[2, i] = 1
    cone_coeffs = VecOrMatOrSparse[coeffs1, coeffs2, coeffs3]
    cone_vars = [x.uid, y.uid, this.uid]
    cone_constant = zeros(3, 1)
    push!(canon_constr_array, CanonicalConstr(cone_coeffs, cone_vars, cone_constant, false, true))
  end

  # y >= 0 linear constraint
  lin_coeffs = VecOrMatOrSparse[-speye(x_size)]
  lin_vars = [y.uid]
  lin_constant = zeros(x_size, 1)
  push!(canon_constr_array, CanonicalConstr(lin_coeffs, lin_vars, lin_constant, false, false))
  append!(canon_constr_array, x.canon_form())
  append!(canon_constr_array, y.canon_form())
  this.canon_form = ()->canon_constr_array
  this.evaluate = ()->x.evaluate() .^ 2 ./ y.evaluate()

  return this
end

qol_elementwise(x::Value, y::Value) = qol_elementwise(convert(CvxExpr, x), convert(CvxExpr, y))
qol_elementwise(x::Value, y::AbstractCvxExpr) = qol_elementwise(convert(CvxExpr, x), y)
qol_elementwise(x::AbstractCvxExpr, y::Value) = qol_elementwise(x, convert(CvxExpr, y))
