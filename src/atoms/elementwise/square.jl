import Base.square
export square

# By using qol_elementwise, square is implemented as x.^2 ./ ones(...)
function square(x::AbstractCvxExpr)
  return qol_elementwise(x, ones(x.size...))
end

# TODO: Extend to generic powers
function ^(x::AbstractCvxExpr, n::Integer)
  if n != 2
    error("Only x^2 is supported as of now")
  end
  return square(x)
end

.^(x::AbstractCvxExpr, n::Integer) = ^(x, n)
