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
    Convex.conv(x::Convex.AbstractExpr, y::Convex.AbstractExpr)

The convolution between two vectors `x` and `y`.

## Examples

```jldoctest; filter=r"id: [0-9]+…[0-9]+"
julia> x = Variable(2);

julia> y = [2, 4];

julia> atom = conv(x, y)
* (affine; real)
├─ 3×2 SparseArrays.SparseMatrixCSC{Int64, Int64} with 4 stored entries
└─ 2-element real variable (id: 663…363)

julia> size(atom)
(3, 1)
```
"""
function conv(x::Value, y::AbstractExpr)
    if length(x) != size(x, 1) || size(y, 2) > 1
        error("convolution only supported between two vectors")
    elseif length(x) == 0
        throw(ArgumentError("convolution with empty vector not supported"))
    end
    X = conv1D_matrix(x, length(y))
    return X * y
end

conv(x::AbstractExpr, y::Value) = conv(y, x)
