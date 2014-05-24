export neg

function neg(x::AbstractCvxExpr)
  this = max(-x, 0)
  this.sign = :neg
  return this
end
