export quadform

function quadform(x::Value, A::AbstractExpr)
  return x' * A * x
end

function quadform(x::AbstractExpr, A::Value)
  if length(size(A)) != 2 || size(A, 1) != size(A, 2)
    error("Quadratic form only takes square matrices")
  end
  if !issymmetric(A)
    error("Quadratic form only defined for symmetric matrices")
  end
  V = eigvals(Symmetric(full(A)))

  if all(V .>= 0)
    factor = 1
  elseif all(V .<= 0)
    factor = -1
  else
    error("Quadratic forms supported only for semidefinite matrices")
  end

  P = real(sqrt(full(factor * A)))
  return factor * square(norm2(P * x))
end
