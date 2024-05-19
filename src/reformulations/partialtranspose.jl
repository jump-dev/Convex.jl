# Copyright (c) 2014: Madeleine Udell and contributors
#
# Use of this source code is governed by a BSD-style license that can be found
# in the LICENSE file or at https://opensource.org/license/bsd-2-clause

"""
    partialtranspose(x, sys::Int, dims::Vector)

Returns the partial transpose of `x` over the `sys`th system, where `dims` is a
vector of integers encoding the dimensions of each subsystem.
"""
function partialtranspose(
    x::Union{AbstractMatrix,AbstractExpr},
    sys::Int,
    dims::Vector,
)
    if size(x, 1) != size(x, 2)
        throw(ArgumentError("Only square matrices are supported"))
    end
    if !(1 <= sys <= length(dims))
        msg = "Invalid system, should between 1 and $(length(dims)); got $sys"
        throw(ArgumentError(msg))
    end
    if size(x, 1) != prod(dims)
        msg = "Dimension of system doesn't correspond to dimension of subsystems"
        throw(ArgumentError(msg))
    end
    n = length(dims)
    s = n - sys + 1
    p = collect(1:2n)
    p[s] = n + s
    p[n+s] = s
    rdims = reverse(dims)
    return _reshape(x, rdims, p)
end

function _reshape(x::AbstractMatrix, rdims, p)
    r = reshape(x, (rdims..., rdims...))
    return reshape(permutedims(r, p), size(x))
end

function _reshape(x::AbstractExpr, rdims, p)
    P = permutedims_matrix((rdims..., rdims...), p)
    return reshape(P * vec(x), size(x)...)
end

"""
    permutedims_matrix(dims, p)

Returns a matrix `M` so that for any vector `v` of length `prod(dims)`,
`M*v == vec(permutedims(reshape(v, dims), p))`.
"""
function permutedims_matrix(dims, p)
    d = prod(dims)
    # Generalization of https://stackoverflow.com/a/60680132
    rows = 1:d
    cols = vec(permutedims(reshape(rows, dims), p))
    data = ones(Int, d)
    return SparseArrays.sparse(rows, cols, data, d, d)
end
