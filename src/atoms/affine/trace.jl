import LinearAlgebra.tr
export trace, tr

function trace(e::AbstractExpr)
  return sum(diag(e))
end

tr(e::AbstractExpr) = trace(e)
