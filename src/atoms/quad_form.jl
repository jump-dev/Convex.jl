export quad_form

function quad_form(x::Constant, A::Constant)
  return x' * A * x
end

function quad_form(x::Constant, A::AbstractCvxExpr)
  return x' * A * x
end

# Note that we don't use Base.issym due to accurracy issues in issym
function is_symmetric(A::Matrix)
  TOLERANCE = 1e-6
  return all(A - A' .< TOLERANCE)
end

function quad_form(x::AbstractCvxExpr, A::Constant)
  if A.size[1] != A.size[2]
    error("Quadratic form only takes square matrices")
  end
  if !is_symmetric(A.value)
    error("Quadratic form only defined for symmetric matrices")
  end
  V = eigvals(full(A.value))
  if !all(V .>= 0) && !all(V .<= 0)
    error("Quadratic forms supported only for semidefinite matrices")
  end

  if all(V .>= 0)
    factor = 1
  else
    factor = -1
  end

  P = sqrtm(full(factor*A.value))
  return factor * square(norm_2(P * x))
end

quad_form(x::Value, A::Value) = quad_form(convert(CvxExpr, x), convert(CvxExpr, A))
quad_form(x::Value, A::AbstractCvxExpr) = quad_form(convert(CvxExpr, x), A)
quad_form(x::AbstractCvxExpr, A::Value) = quad_form(x, convert(CvxExpr, A))
