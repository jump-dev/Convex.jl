# Copyright (c) 2014: Madeleine Udell and contributors
#
# Use of this source code is governed by a BSD-style license that can be found
# in the LICENSE file or at https://opensource.org/license/bsd-2-clause

"""
    Base.kron(x::Convex.AbstractExpr, y::Convex.AbstractExpr)

The Kronecker (outer) product.

## Examples

```jldoctest; filter=r"id: [0-9]+…[0-9]+"
julia> x = Variable(2);

julia> y = [1 2];

julia> atom = kron(x, y)
vcat (affine; real)
├─ * (affine; real)
│  ├─ index (affine; real)
│  │  └─ 2-element real variable (id: 369…232)
│  └─ [1 2]
└─ * (affine; real)
   ├─ index (affine; real)
   │  └─ 2-element real variable (id: 369…232)
   └─ [1 2]

julia> size(atom)
(2, 2)
```
"""
function Base.kron(a::Value, b::AbstractExpr)
    rows = AbstractExpr[]
    for i in 1:size(a)[1]
        row = AbstractExpr[]
        for j in 1:size(a)[2]
            if isreal(a[i, j])
                push!(row, real(a[i, j]) * b)
            else
                push!(row, a[i, j] * b)
            end
        end
        push!(rows, foldl(hcat, row))
    end
    return foldl(vcat, rows)
end

function Base.kron(a::AbstractExpr, b::Value)
    rows = AbstractExpr[]
    for i in 1:size(a)[1]
        row = AbstractExpr[]
        for j in 1:size(a)[2]
            if isreal(b)
                push!(row, a[i, j] * real(b))
            else
                push!(row, a[i, j] * b)
            end
        end
        push!(rows, foldl(hcat, row))
    end
    return foldl(vcat, rows)
end
