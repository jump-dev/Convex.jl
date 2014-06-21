export pos

function pos(x::AbstractCvxExpr)
  this = max(x, 0)
  this.sign = :pos
  return this
end
