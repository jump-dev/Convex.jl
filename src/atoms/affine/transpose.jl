function LinearAlgebra.transpose(x::AbstractExpr)
    sz = length(x)

    num_rows = x.size[1]
    num_cols = x.size[2]

    I = Array{Int}(undef, sz)
    J = Array{Int}(undef, sz)

    k = 1
    for r in 1:num_rows
        for c in 1:num_cols
            I[k] = (c - 1) * num_rows + r
            J[k] = (r - 1) * num_cols + c
            k += 1
        end
    end

    transpose_matrix = sparse(I, J, 1.0)

    return reshape(transpose_matrix * vec(x), size(x, 2), size(x, 1))
end

function LinearAlgebra.transpose(x::Union{Constant,ComplexConstant})
    return constant(transpose(x.value))
end

LinearAlgebra.adjoint(x::AbstractExpr) = transpose(conj(x))
LinearAlgebra.adjoint(x::Union{Constant,ComplexConstant}) = constant(x.value')
