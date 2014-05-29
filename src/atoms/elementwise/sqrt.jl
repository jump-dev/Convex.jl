import Base.sqrt
export sqrt

function sqrt(x::AbstractCvxExpr)
  return geo_mean(x, ones(x.size...))
end
