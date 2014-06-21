export qol_elementwise

# The qol_elementwise atom takes two expressions, x and y, of the same size
# and returns an expression of the same size, where each element is given by
# x_i^2 / y_i. Additionally, it is required that y_i >= 0 for all i.

function qol_elementwise(x::Constant, y::Constant)
  return Constant((x.value.^2)./y.value, :pos)
end

function qol_elementwise(x::Constant, y::AbstractCvxExpr)
  if x.size != y.size
    error("Expressions must have same size for qol_elementwise")
  end
  if !is_concave(y.vexity)
    error("Denominator of qol_elementwise cannot be a convex expression")
  end
  this = CvxExpr(:qol_elementwise, [x, y], :convex, :pos, x.size)
  x_size = get_vectorized_size(x)

  # If t is the new expression created, then the canonical form is given by
  # the SOCP constraint (y_i + t_i, y_i - t_i, 2 * x_i) \in Q_3
  # for each index i. This is equivalent to the conic inequality
  #   0  0     y_i      2 * x_i
  #  -1 -1  *  t_i  <=  0
  #  -1  1              0
  # which is constructed by the following code
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

  # Additionally, y_i >= 0 is also enforced
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
  if x.size != y.size
    error("Expressions must have same size for qol_elementwise")
  end
  if x.sign == :pos && !is_convex(x.vexity) ||
     x.sign == :neg && !is_concave(x.vexity) ||
     x.sign == :any && !is_linear(x.vexity)
    error("Numerator of qol_elementwise expression is not DCP compliant")
  elseif y.sign == :neg
    error("Denominator of qol_elementwise cannot be a negative constant")
  end

  this = CvxExpr(:qol_elementwise, [x, y], :convex, :pos, x.size)
  x_size = get_vectorized_size(x)

  # If t is the new expression created, then the canonical form is given by
  # the SOCP constraint (y_i + t_i, y_i - t_i, 2 * x_i) \in Q_3
  # for each index i. This is equivalent to the conic inequality
  #  -2  0     x_i      0
  #   0 -1  *  t_i  <=  y_i
  #   0  1              y_i
  # which is constructed by the following code
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
  if x.size != y.size
    error("Expressions must have same size for qol_elementwise")
  end
  if x.sign == :pos && !is_convex(x.vexity) ||
     x.sign == :neg && !is_concave(x.vexity) ||
     x.sign == :any && !is_affine(x.vexity)
    error("Numerator of qol_elementwise expression is not DCP compliant")
  elseif !is_concave(y.vexity)
    error("Denominator of qol_elementwise cannot be a convex expression")
  end
  this = CvxExpr(:qol_elementwise, [x, y], :convex, :pos, x.size)
  x_size = get_vectorized_size(x)

  # If t is the new expression created, then the canonical form is given by
  # the SOCP constraint (y_i + t_i, y_i - t_i, 2 * x_i) \in Q_3
  # for each index i. This is equivalent to the conic inequality
  #  -2  0  0     x_i      0
  #   0 -1 -1  *  y_i  <=  0
  #   0 -1  1     t_i      0
  # which is constructed by the following code
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

  # Additionally, y_i >= 0 is also enforced
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
