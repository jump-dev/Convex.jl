# Copyright (c) 2014: Madeleine Udell and contributors
#
# Use of this source code is governed by a BSD-style license that can be found
# in the LICENSE file or at https://opensource.org/license/bsd-2-clause

_eps_tol(A) = sqrt(eps(float(real(eltype(A)))))

"""
    _is_psd(A; tol)

Check whether `A` is positive semi-definite by computing a LDLᵀ factorization of
`A + tol*I`
"""
function _is_psd(A; tol = _eps_tol(A))
    T = eltype(A)
    # If `A` is neither a Matrix nor SparseMatrixCSC, we do the following:
    # * sparse fallback if the arithmetic is supported
    # * dense fallack otherwise
    if T <: AbstractFloat || T <: LinearAlgebra.BlasFloat
        return _is_psd(SparseArrays.sparse(A); tol = tol)
    end
    return _is_psd(Matrix(A); tol = tol)
end

function _is_psd(A::Matrix; tol = _eps_tol(A))
    return LinearAlgebra.isposdef(A + tol * LinearAlgebra.I)
end

function _is_psd(
    A::SparseArrays.SparseMatrixCSC{Complex{T}};
    tol::T = _eps_tol(A),
) where {T<:LinearAlgebra.BlasReal}
    return LinearAlgebra.isposdef(A + tol * LinearAlgebra.I)
end

function _is_psd(
    A::SparseArrays.SparseMatrixCSC{T};
    tol::T = _eps_tol(A),
) where {T<:AbstractFloat}
    # LDLFactorizations requires the input matrix to only have the upper triangle.
    A_ = LinearAlgebra.Symmetric(
        SparseArrays.sparse(LinearAlgebra.UpperTriangular(A)) +
        tol * LinearAlgebra.I,
    )
    try
        F = LDLFactorizations.ldl(A_)
        return minimum(F.D.diag) >= 0
    catch err
        # If the matrix could not be factorized, then it is not PSD
        if err isa LDLFactorizations.SQDException
            return false
        end
        rethrow(err)
    end
end

"""
    quadform(x::AbstractExpr, A::AbstractExpr; assume_psd=false)

Represents \$x^\\top A x\$ where either:

 * `x` is a vector-valued variable and `A` is a positive semidefinite or
   negative semidefinite matrix (and in particular Hermitian or real symmetric).
   If `assume_psd=true`, then `A` will be assumed to be positive semidefinite.
   Otherwise, `Convex._is_psd` will be used to check if `A` is positive
   semidefinite or negative semidefinite.
 * or `A` is a matrix-valued variable and `x` is a vector.

## Examples

```jldoctest; filter=r"id: [0-9]+…[0-9]+"
julia> x = Variable(2);

julia> A = [1 0; 0 1]
2×2 Matrix{Int64}:
 1  0
 0  1

julia> atom = quadform(x, A)
* (convex; positive)
├─ [1;;]
└─ qol_elem (convex; positive)
   ├─ norm2 (convex; positive)
   │  └─ * (affine; real)
   │     ├─ …
   │     └─ …
   └─ [1.0;;]

julia> size(atom)
(1, 1)
```

```jldoctest; filter=r"id: [0-9]+…[0-9]+"
julia> x = [1, 2]

julia> A = Variable(2, 2);

julia> atom = quadform(x, A)
* (affine; real)
├─ * (affine; real)
│  ├─ [1 2]
│  └─ 2×2 real variable (id: 111…794)
└─ [1; 2;;]

julia> size(atom)
(1, 1)
```
"""
quadform(x::Value, A::AbstractExpr; kwargs...) = x' * A * x

function quadform(x::AbstractExpr, A::Value; assume_psd = false)
    if length(size(A)) != 2 || size(A, 1) != size(A, 2)
        error("quadform only takes square matrices")
    elseif !LinearAlgebra.ishermitian(A)
        error("quadform only defined for Hermitian matrices")
    end
    factor = if assume_psd || _is_psd(A)
        1
    elseif _is_psd(-A)
        -1
    else
        error("Quadratic forms supported only for semidefinite matrices")
    end
    # The term `factor * A` is SPSD
    P = _squareroot(factor * A)
    return factor * square(norm2(P * x))
end

function _squareroot(x::AbstractExpr, A::Value;)
    # Caluclates matrix `P` such that `A = P' * P`
    P = sqrt(LinearAlgebra.Hermitian(A))
    return P
end

function _squareroot(x::AbstractExpr, A::SparseArrays.SparseMatrixCSC; shift = 1e-10)
    # Caluclates matrix `P` such that `A = P' * P`
    _chol = cholesky(A; shift = shift, check = false)
    if !issuccess(_chol)
        error("The Shifted Cholesky decomposition failed. The matrix may not be an SPSD.");
    end
    # Extract sparse L and permutation
    L = sparse(_chol.L);
    p = _chol.p;
    pinv = invperm(p);
    P = L'[:, pinv];
    return P
end

# These `evaluate` calls are safe since they are not a `fix!`ed variable
# (must be a constant):
function quadform(x::Constant, A::AbstractExpr; kwargs...)
    return quadform(evaluate(x), A; kwargs...)
end

function quadform(x::AbstractExpr, A::Constant; kwargs...)
    return quadform(x, evaluate(A); kwargs...)
end

function quadform(x::Constant, A::Constant; kwargs...)
    return evaluate(x)' * evaluate(A) * evaluate(x)
end

function quadform(x::AbstractExpr, A::AbstractExpr; kwargs...)
    return error(
        "Convex.jl v0.13.5 introduced the ability to use `fix!`ed variables " *
        "in `quadform`. However, this did not consider the case that the " *
        "value of `fix!`ed variables is changed between solves. Due to the " *
        "risk that this may silently produce incorrect solutions, this " *
        "behavior has been removed. Use `evaluate(H)` to obtain the value of " *
        "a fixed variable. If the value changes between solves, rebuild the " *
        "problem for the change to take effect.",
    )
end
