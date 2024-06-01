# Copyright (c) 2014: Madeleine Udell and contributors
#
# Use of this source code is governed by a BSD-style license that can be found
# in the LICENSE file or at https://opensource.org/license/bsd-2-clause

"""
    conv1D_matrix(h::AbstractVector, n::Integer) -> SparseMatrixCSC

Create a sparse matrix `A` such that if `x` has length `n`,
then we have `A * x ≈ conv1d(h, x)`.
"""
function conv1D_matrix(h::AbstractVector, n::Integer)
    m = length(h)
    Is = Int[]
    Js = Int[]
    Vs = eltype(h)[]
    sizehint!(Is, n * m)
    sizehint!(Js, n * m)
    sizehint!(Vs, n * m)
    # build matrix by columns
    for j in 1:n
        append!(Is, j:(j+m-1))
        append!(Js, (j for _ in 1:m))
        append!(Vs, h)
    end
    return SparseArrays.sparse(Is, Js, Vs, m + n - 1, n)
end

"""
    conv(x, y)

Produces the discrete convolution of vectors `x` of length `m` and `y` of length `n`. That is, if `z = conv(x,y)`, then `z` has length `m+n-1`, and `z` has entries

```julia
z[i] = sum(x[j] * get(y, i - j + 1, 0) for j in 1:m)
```

Note that `conv` is symmetric in `x` and `y`: `conv(x,y) == conv(y, x)` (in exact arithmetic).
"""
function conv(x::Value, y::AbstractExpr)
    if length(x) != size(x, 1) || size(y, 2) > 1
        error("convolution only supported between two vectors")
    end
    if length(x) == 0
        throw(ArgumentError("convolution with empty vector not supported"))
    end
    X = conv1D_matrix(x, length(y))
    return X * y
end

conv(x::AbstractExpr, y::Value) = conv(y, x)

# direct non-variable implementation for reference
function conv(x::AbstractVector, y::AbstractVector)
    T = promote_type(eltype(x), eltype(y))
    m = length(x)
    n = length(y)
    z = zeros(T, m + n - 1)
    for i in eachindex(z)
        z[i] = sum(x[j] * get(y, i - j + 1, 0) for j in 1:m)
    end
    return z
end

#####
##### 2D discrete convolution
#####

# We reformulate the problem into a 1D convolution following the approach from:
# https://en.wikipedia.org/wiki/Multidimensional_discrete_convolution#Multidimensional_convolution_with_one-dimensional_convolution_methods

# In particular, we keep separate the matrix representation of the 1D convolution, since this could be re-used among multiple 2D convolutions.

function conv2D_kernel(X::AbstractMatrix, sz::Tuple{Int,Int})
    M, N = size(X)
    K, L = sz
    Z_rows = M + K - 1
    Z_cols = N + L - 1
    Xpad = zeros(eltype(X), Z_rows, Z_cols)
    Xpad[1:size(X, 1), 1:size(X, 2)] .= X
    Xv = vec(Xpad)[1:(N-1)*(M+K-1)+M]
    return conv1D_matrix(Xv, (L - 1) * (M + K - 1) + K)
end

function conv2D(
    X_size::Tuple,
    X_kernel::AbstractMatrix{T},
    Y::Convex.AbstractExpr,
) where {T}
    M, N = X_size
    K, L = size(Y)
    Z_rows = M + K - 1
    Z_cols = N + L - 1
    bottom = spzeros(T, length(size(Y, 1)+1:Z_rows), Z_cols)
    right_side = spzeros(T, size(Y, 1), length(size(Y, 2)+1:Z_cols))
    Y_padded = [
        Y right_side
        bottom
    ]
    Y_vector = vec(Y_padded)[1:(L-1)*(M+K-1)+K]
    Z_vector = X_kernel * Y_vector
    Z = reshape(Z_vector, Z_rows, Z_cols)
    return Z
end

"""
    conv2D(X, Y)

Performs a 2D discrete convolution of `X` and `Y`. Assuming `size(X) = (M,N)` and `size(Y) = (K, L)`, then `Z = conv2D(X,Y)` has `size(Z) = (M + K - 1, N + L - 1)`, with entries

```julia
Z[i, j] = sum(X[m, n] * get(Y, (i-m+1, j-n+1), zero(T)) for m in 1:M, n in 1:N)
```

"""
function conv2D(X::AbstractMatrix, Y::Convex.AbstractExpr)
    return conv2D(size(X), conv2D_kernel(X, size(Y)), Y)
end

conv2D(Y::Convex.AbstractExpr, X::AbstractMatrix) = conv2D(X, Y)

# direct non-variable implementation for reference
function conv2D(X::AbstractMatrix, Y::AbstractMatrix)
    T = promote_type(eltype(X), eltype(Y))
    M, N = size(X)
    K, L = size(Y)
    Z = zeros(T, M + K - 1, N + L - 1)
    for i in 1:M+K-1, j in 1:N+L-1
        Z[i, j] = sum(
            X[m, n] * get(Y, (i - m + 1, j - n + 1), zero(T)) for m in 1:M,
            n in 1:N
        )
    end
    return Z
end
