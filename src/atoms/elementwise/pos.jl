export pos

function pos(x::AbstractCvxExpr)
  this = max(x, spzeros(x.size...))
  this.sign = :pos
  return this
end
