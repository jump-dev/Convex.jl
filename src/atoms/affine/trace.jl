import LinearAlgebra.tr
export tr

function tr(e::AbstractExpr)
    return sum(diag(e))
end
