# Since everything is vectorized, we simply need to multiply x by a permutation
# matrix such that coeff * vectorized(x) - vectorized(transpose(x)) = 0
function LinearAlgebra.transpose(x::AbstractExpr)
    transpose_matrix = permutedims_matrix(size(x), (2, 1))
    return reshape(transpose_matrix * vec(x), size(x, 2), size(x, 1))
end

function LinearAlgebra.transpose(x::Union{Constant,ComplexConstant})
    return constant(transpose(x.value))
end

LinearAlgebra.adjoint(x::AbstractExpr) = transpose(conj(x))
LinearAlgebra.adjoint(x::Union{Constant,ComplexConstant}) = constant(x.value')
