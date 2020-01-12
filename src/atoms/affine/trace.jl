import LinearAlgebra.tr

function tr(e::AbstractExpr)
    return sum(diag(e))
end
