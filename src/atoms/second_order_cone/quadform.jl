function quadform(x::Value, A::AbstractExpr; assume_psd=false)
    return x' * A * x
end

"""
    is_psd(A; tol)

Check whether `A` is positive semi-definite by computing a LDLᵀ factorization of `A + tol*I`
"""
function is_psd(A; tol=sqrt(eps(float(real(eltype(A))))))
    T = eltype(A)
    # If `A` is neither a Matrix nor SparseMatrixCSC, we do the following:
    # * sparse fallback if the arithmetic is supported
    # * dense fallack otherwise
    if T <: Real || T <: LinearAlgebra.BlasFloat
        return is_psd(sparse(A); tol=tol)
    else
        return is_psd(Matrix(A); tol=tol)
    end
end
is_psd(A::SparseMatrixCSC{Complex{T}}; tol::T=sqrt(eps(T))) where{T<:LinearAlgebra.BlasReal} = isposdef(A + tol*I)
function is_psd(A::SparseMatrixCSC{T}; tol::T=sqrt(eps(T))) where{T<:Real}
    # LDLFactorizations requires the input matrix to only have the upper triangle.
    A_ = Symmetric(sparse(UpperTriangular(A)) + tol*I)
    try
        F = ldl(A_)
        d = F.D.diag
        return minimum(d) >= 0
    catch err
        # If the matrix could not be factorized, then it is not PSD
        isa(err, LDLFactorizations.SQDException) && return false
        rethrow()  # Something else happened
    end
end
is_psd(A::Matrix; tol=sqrt(eps(float(real(eltype(A)))))) = isposdef(A + tol*I)

function quadform(x::AbstractExpr, A::Value; assume_psd=false)
    if length(size(A)) != 2 || size(A, 1) != size(A, 2)
        error("Quadratic form only takes square matrices")
    end
    if !ishermitian(A)
        error("Quadratic form only defined for Hermitian matrices")
    end
    if assume_psd
        factor = 1
    else
        if is_psd(A)
            factor = 1
        elseif is_psd(-A)
            factor = -1
        else
            error("Quadratic forms supported only for semidefinite matrices")
        end
    end

    P = sqrt(Hermitian(factor * A))
    return factor * square(norm2(P * x))
end

"""
    quadform(x::AbstractExpr, A::AbstractExpr; assume_psd=false)

Represents `x' * A * x` where either:

* `x` is a vector-valued variable and `A` is a positive semidefinite or negative semidefinite matrix (and in particular Hermitian or real symmetric). If `assume_psd=true`, then
`A` will be assumed to be positive semidefinite. Otherwise, `Convex.is_psd` will be used to check if `A` is positive semidefinite or negative semidefinite.
* or `A` is a matrix-valued variable and `x` is a vector.

"""
function quadform(x::AbstractExpr, A::AbstractExpr; assume_psd=false)
    if vexity(x) == ConstVexity()
        return quadform(evaluate(x), A; assume_psd=assume_psd)
    elseif vexity(A) == ConstVexity()
        return quadform(x, evaluate(A); assume_psd=assume_psd)
    else
        error("Either `x` or `A` must be constant in `quadform(x,A)`.")
    end
end
