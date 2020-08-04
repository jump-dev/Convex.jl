function LinearAlgebra.tr(e::AbstractExpr)
    return sum(diag(e))
end
