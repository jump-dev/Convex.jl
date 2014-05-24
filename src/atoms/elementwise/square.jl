export square

function square(x::AbstractCvxExpr)
  return qol_elementwise(x, ones(x.size...))
end
