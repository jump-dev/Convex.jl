function quadform(x::Value, A::AbstractExpr; assume_psd = false)
    return x' * A * x
end

"""
    is_psd(A; tol)

Check whether `A` is positive semi-definite by computing a LDLᵀ factorization of
`A + tol*I`
"""
function is_psd(A; tol = sqrt(eps(float(real(eltype(A))))))
    T = eltype(A)
    # If `A` is neither a Matrix nor SparseMatrixCSC, we do the following:
    # * sparse fallback if the arithmetic is supported
    # * dense fallack otherwise
    if T <: AbstractFloat || T <: LinearAlgebra.BlasFloat
        return is_psd(SparseArrays.sparse(A); tol = tol)
    end
    return is_psd(Matrix(A); tol = tol)
end

function is_psd(
    A::SparseArrays.SparseMatrixCSC{Complex{T}};
    tol::T = sqrt(eps(T)),
) where {T<:LinearAlgebra.BlasReal}
    return LinearAlgebra.isposdef(A + tol * LinearAlgebra.I)
end

function is_psd(
    A::SparseArrays.SparseMatrixCSC{T};
    tol::T = sqrt(eps(T)),
) where {T<:AbstractFloat}
    # LDLFactorizations requires the input matrix to only have the upper triangle.
    A_ = LinearAlgebra.Symmetric(
        SparseArrays.sparse(LinearAlgebra.UpperTriangular(A)) +
        tol * LinearAlgebra.I,
    )
    try
        F = LDLFactorizations.ldl(A_)
        d = F.D.diag
        return minimum(d) >= 0
    catch err
        # If the matrix could not be factorized, then it is not PSD
        if err isa LDLFactorizations.SQDException
            return false
        end
        rethrow(err)
    end
end

function is_psd(A::Matrix; tol = sqrt(eps(float(real(eltype(A))))))
    return LinearAlgebra.isposdef(A + tol * LinearAlgebra.I)
end

function quadform(x::AbstractExpr, A::Value; assume_psd = false)
    if length(size(A)) != 2 || size(A, 1) != size(A, 2)
        error("Quadratic form only takes square matrices")
    elseif !LinearAlgebra.ishermitian(A)
        error("Quadratic form only defined for Hermitian matrices")
    end
    factor = if assume_psd || is_psd(A)
        1
    elseif is_psd(-A)
        -1
    else
        error("Quadratic forms supported only for semidefinite matrices")
    end
    P = sqrt(LinearAlgebra.Hermitian(factor * A))
    return factor * square(norm2(P * x))
end

"""
    quadform(x::AbstractExpr, A::AbstractExpr; assume_psd=false)

Represents `x' * A * x` where either:

 * `x` is a vector-valued variable and `A` is a positive semidefinite or
   negative semidefinite matrix (and in particular Hermitian or real symmetric).
   If `assume_psd=true`, then `A` will be assumed to be positive semidefinite.
   Otherwise, `Convex.is_psd` will be used to check if `A` is positive
   semidefinite or negative semidefinite.
 * or `A` is a matrix-valued variable and `x` is a vector.
"""
function quadform(x::AbstractExpr, A::AbstractExpr; assume_psd = false)
    if vexity(x) == ConstVexity()
        return quadform(evaluate(x), A; assume_psd = assume_psd)
    elseif vexity(A) == ConstVexity()
        return quadform(x, evaluate(A); assume_psd = assume_psd)
    else
        error("Either `x` or `A` must be constant in `quadform(x,A)`.")
    end
end