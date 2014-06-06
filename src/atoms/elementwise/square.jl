import Base.square
export square

# By using qol_elementwise, square is implemented as x.^2 ./ ones(...)
function square(x::AbstractCvxExpr)
  return qol_elementwise(x, ones(x.size...))
end
