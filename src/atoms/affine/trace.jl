import Base.trace
export trace

function trace(e::AbstractExpr)
  return sum(diag(e))
end
