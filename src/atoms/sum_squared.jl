export sum_squared

function sum_squared(x::AbstractCvxExpr)
  return square(norm_2(x))
end
