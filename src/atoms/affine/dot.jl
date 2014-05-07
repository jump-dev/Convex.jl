export dot

function dot(x::AbstractCvxExpr, y::AbstractCvxExpr)
  return x' * y
end

dot(x::Value, y::AbstractCvxExpr) = dot(convert(CvxExpr, x), y)
dot(x::AbstractCvxExpr, y::Value) = dot(x, convert(CvxExpr, y))
