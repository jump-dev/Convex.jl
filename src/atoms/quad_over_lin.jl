export quad_over_lin

# TODO: @david this function isn't being used anywhere
function check_size_qol(x::AbstractCvxExpr, y::AbstractCvxExpr)
  if (x.size[1] > 1 && x.size[2] > 1) || y.size != (1, 1)
    error("quad_over_lin arguments must be a vector and a scalar")
  end
end

function quad_over_lin(x::Constant, y::Constant)
  #TODO sign/size checks
  return x' * x / y
end

function quad_over_lin(x::Constant, y::AbstractCvxExpr)
  #TODO vexity and sign checks
  this = CvxExpr(:quad_over_lin, [x, y], :convex, :pos, (1, 1))
  x_size = get_vectorized_size(x)
  cone_size = x_size + 2

  # (y + t, y - t, 2x) socp constraint
  coeffs1 = spzeros(cone_size, 1)
  coeffs1[1] = -1
  coeffs1[2] = -1
  coeffs2 = spzeros(cone_size, 1)
  coeffs2[1] = -1
  coeffs2[2] = 1
  cone_coeffs = VecOrMatOrSparse[coeffs1, coeffs2]
  cone_vars = [y.uid, this.uid]
  cone_constant = [0; 0; 2*vec(x.value)]

  # y >= 0 linear constraint
  lin_coeffs = VecOrMatOrSparse[-speye(1)]
  lin_vars = [y.uid]
  lin_constant = zeros(1, 1)

  canon_constr_array = [CanonicalConstr(cone_coeffs, cone_vars, cone_constant, false, true),
                        CanonicalConstr(lin_coeffs, lin_vars, lin_constant, false, false)]
  append!(canon_constr_array, y.canon_form())
  this.canon_form = ()->canon_constr_array
  this.evaluate = ()->x.evaluate()' * x.evaluate() / y.evaluate()

  return this
end

function quad_over_lin(x::AbstractCvxExpr, y::Constant)
  #TODO vexity and sign checks
  this = CvxExpr(:quad_over_lin, [x, y], :convex, :pos, (1, 1))
  x_size = get_vectorized_size(x)
  cone_size = x_size + 2

  # (y + t, y - t, 2x) socp constraint
  coeffs1 = [spzeros(2, x_size); -2*speye(x_size)]
  coeffs2 = spzeros(cone_size, 1)
  coeffs2[1] = -1
  coeffs2[2] = 1
  cone_coeffs = VecOrMatOrSparse[coeffs1, coeffs2]
  cone_vars = [x.uid, this.uid]
  cone_constant = [y.value[1]; y.value[1]; zeros(x_size, 1)]

  canon_constr_array = [CanonicalConstr(cone_coeffs, cone_vars, cone_constant, false, true)]
  append!(canon_constr_array, x.canon_form())
  this.canon_form = ()->canon_constr_array
  this.evaluate = ()->x.evaluate()' * x.evaluate() / y.evaluate()

  return this
end

function quad_over_lin(x::AbstractCvxExpr, y::AbstractCvxExpr)
  #TODO vexity and sign checks
  this = CvxExpr(:quad_over_lin, [x, y], :convex, :pos, (1, 1))
  x_size = get_vectorized_size(x)
  cone_size = x_size + 2

  # (y + t, y - t, 2x) socp constraint
  coeffs1 = [spzeros(2, x_size); -2*speye(x_size)]
  coeffs2 = spzeros(cone_size, 1)
  coeffs2[1] = -1
  coeffs2[2] = -1
  coeffs3 = spzeros(cone_size, 1)
  coeffs3[1] = -1
  coeffs3[2] = 1
  cone_coeffs = VecOrMatOrSparse[coeffs1, coeffs2, coeffs3]
  cone_vars = [x.uid, y.uid, this.uid]
  cone_constant = zeros(cone_size, 1)

  # y >= 0 linear constraint
  lin_coeffs = VecOrMatOrSparse[-speye(1)]
  lin_vars = [y.uid]
  lin_constant = zeros(1, 1)

  canon_constr_array = [CanonicalConstr(cone_coeffs, cone_vars, cone_constant, false, true),
                        CanonicalConstr(lin_coeffs, lin_vars, lin_constant, false, false)]
  append!(canon_constr_array, x.canon_form())
  append!(canon_constr_array, y.canon_form())
  this.canon_form = ()->canon_constr_array

  this.evaluate = ()->x.evaluate()' * x.evaluate() / y.evaluate()

  return this
end

quad_over_lin(x::Value, y::Value) = quad_over_lin(convert(CvxExpr, x), convert(CvxExpr, y))
quad_over_lin(x::Value, y::AbstractCvxExpr) = quad_over_lin(convert(CvxExpr, x), y)
quad_over_lin(x::AbstractCvxExpr, y::Value) = quad_over_lin(x, convert(CvxExpr, y))
