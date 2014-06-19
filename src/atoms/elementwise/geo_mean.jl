export geo_mean

# Find the elementwise geometric mean between two expressions x and y
# e.g. for two vectors x and y, compute
# (sqrt(x_1 * y_1), sqrt(x_2 * y_2), ..., sqrt(x_n * y_n))
# TODO: We don't allow geometric mean along a dimension of an expression yet

function geo_mean(x::Constant, y::Constant)
  return Constant(sqrt(x.value .* y.value), :pos)
end

function geo_mean(x::AbstractCvxExpr, y::Constant)
  if x.size != y.size
    error("Geometric mean only supported between expressions of the same size")
  end
  if !is_concave(x.vexity)
    error("Geometric mean of convex expressions is not DCP compliant")
  end
  if y.sign != :pos
    error("Cannot take geometric mean of non-positive constant")
  end
  this = CvxExpr(:geo_mean, [x, y], :concave, :pos, x.size)
  x_size = get_vectorized_size(x)

  # If t is the new expression created, then the canonical form is given by
  # the SOCP constraint (y_i + x_i, y_i - x_i, 2 * t_i) \in Q_3
  # for each index i. This is equivalent to the conic inequality
  #  -1  0     x_i      y_i
  #   1  0  *  t_i  <=  y_i
  #   0 -2              0
  # which is constructed by the following code
  canon_constr_array = CanonicalConstr[]
  for i = 1:x_size
    coeffs1 = spzeros(3, x_size)
    coeffs1[1, i] = -1
    coeffs1[2, i] = 1
    coeffs2 = spzeros(3, x_size)
    coeffs2[3, i] = -2
    cone_coeffs = VecOrMatOrSparse[coeffs1, coeffs2]
    cone_vars = [x.uid, this.uid]
    cone_constant = [y.value[i], y.value[i], 0]
    push!(canon_constr_array, CanonicalConstr(cone_coeffs, cone_vars, cone_constant, false, true))
  end

  # In addition, x_i >= 0 must also be enforced
  lin_coeffs = VecOrMatOrSparse[-speye(x_size)]
  lin_vars = [x.uid]
  lin_constant = zeros(x_size, 1)
  push!(canon_constr_array, CanonicalConstr(lin_coeffs, lin_vars, lin_constant, false, false))
  append!(canon_constr_array, x.canon_form())
  this.canon_form = ()->canon_constr_array
  this.evaluate = ()->Base.sqrt(x.evaluate() .* y.evaluate())

  return this
end

geo_mean(x::Constant, y::AbstractCvxExpr) = geo_mean(y, x)

function geo_mean(x::AbstractCvxExpr, y::AbstractCvxExpr)
  if x.size != y.size
    error("Geometric mean only supported between expressions of the same size")
  end
  if !is_concave(x.vexity) || !is_concave(y.vexity)
    errory("Geometric mean of convex expressions is not DCP compliant")
  end
  this = CvxExpr(:geo_mean, [x, y], :concave, :pos, x.size)
  x_size = get_vectorized_size(x)

  # If t is the new expression created, then the canonical form is given by
  # the SOCP constraint (y_i + x_i, y_i - x_i, 2 * t_i) \in Q_3
  # for each index i. This is equivalent to the conic inequality
  #  -1 -1  0     x_i      0
  #   1 -1  0  *  y_i  <=  0
  #   0  0 -2     t_i      0
  # which is constructed by the following code
  canon_constr_array = CanonicalConstr[]
  for i = 1:x_size
    coeffs1 = spzeros(3, x_size)
    coeffs1[1, i] = -1
    coeffs1[2, i] = 1
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

  # In addition, x_i, y_i >= 0 must also be enforced
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

geo_mean(x::Value, y::Value) = geo_mean(convert(CvxExpr, x), convert(CvxExpr, y))
geo_mean(x::Value, y::AbstractCvxExpr) = geo_mean(convert(CvxExpr, x), y)
geo_mean(x::AbstractCvxExpr, y::Value) = geo_mean(x, convert(CvxExpr, y))
