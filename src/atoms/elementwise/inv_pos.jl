export inv_pos

function inv_pos(x::AbstractCvxExpr)
  return qol_elementwise(ones(x.size...), x)
end
