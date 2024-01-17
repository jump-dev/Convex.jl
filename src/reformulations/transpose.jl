# Since everything is vectorized, we simply need to multiply x by a permutation
# matrix such that coeff * vectorized(x) - vectorized(transpose(x)) = 0
function LinearAlgebra.transpose(x::AbstractExpr)
    P = permutedims_matrix(size(x), (2, 1))
    return reshape(P * vec(x), size(x, 2), size(x, 1))
end

LinearAlgebra.adjoint(x::AbstractExpr) = transpose(conj(x))

function LinearAlgebra.transpose(x::Union{Constant,ComplexConstant})
    return constant(transpose(evaluate(x)))
end

function LinearAlgebra.adjoint(x::Union{Constant,ComplexConstant})
    return constant(adjoint(evaluate(x)))
end
