export square_pos

function square_pos(x::AbstractCvxExpr)
  return qol_elementwise(pos(x), ones(x.size...))
end
