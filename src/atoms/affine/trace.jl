function LinearAlgebra.tr(e::AbstractExpr)
    return sum(LinearAlgebra.diag(e))
end
