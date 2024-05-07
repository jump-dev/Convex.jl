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
    # It is much more efficient to construct sparse matrices
    # this way rather than starting from `spzeros` and indexing into it.
    Is = Int[]
    Js = Int[]
    Vs = eltype(h)[]
    # build matrix by columns
    for j in 1:n
        append!(Is, j:(j+m-1))
        append!(Js, repeat([j], m))
        append!(Vs, h)
    end
    return SparseArrays.sparse(Is, Js, Vs, m + n - 1, n)
end

function conv(x::Value, y::AbstractExpr)
    if length(x) != size(x, 1) || size(y, 2) > 1
        error("convolution only supported between two vectors")
    end
    length(x) > 0 ||
        throw(ArgumentError("convolution with empty vector not supported"))

    X = conv1D_matrix(x, length(y))
    return X * y
end

conv(x::AbstractExpr, y::Value) = conv(y, x)
