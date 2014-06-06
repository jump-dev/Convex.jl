export square_pos

# Square of pos(x)
function square_pos(x::AbstractCvxExpr)
  return square(pos(x))
end
