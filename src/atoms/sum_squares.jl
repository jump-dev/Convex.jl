export sum_squares

function sum_squares(x::AbstractCvxExpr)
  return square(norm_2(x))
end
