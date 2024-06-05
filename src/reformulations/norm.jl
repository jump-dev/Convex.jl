# Copyright (c) 2014: Madeleine Udell and contributors
#
# Use of this source code is governed by a BSD-style license that can be found
# in the LICENSE file or at https://opensource.org/license/bsd-2-clause

"""
    norm(x::AbstractExpr, p::Real = 2)

Computes the `p`-norm `‖x‖ₚ = (∑ᵢ |xᵢ|^p)^(1/p)` of a vector expression `x`.

Matrices are vectorized (i.e., `norm(x)` is the same as `norm(vec(x))`.)

The return value depends on the value of `p`. Specialized cases are used for
`p = 1`, `p = 2`, and `p = Inf`.

## Examples

```jldoctest; filter=r"id: [0-9]+…[0-9]+"
julia> x = Variable(2);

julia> atom = norm(x, 1)
sum (convex; positive)
└─ abs (convex; positive)
   └─ 2-element real variable (id: 779…899)

julia> size(atom)
(1, 1)

julia> norm(x, 2)
norm2 (convex; positive)
└─ 2-element real variable (id: 779…899)

julia> norm(x, Inf)
maximum (convex; positive)
└─ abs (convex; positive)
   └─ 2-element real variable (id: 779…899)

julia> norm(x, 3 // 2)
rationalnorm (convex; positive)
└─ 2-element real variable (id: 779…899)
```
"""
function LinearAlgebra.norm(x::AbstractExpr, p::Real = 2)
    if size(x, 2) > 1
        x = vec(x)
    end
    if p == 1
        return sum(abs(x))
    elseif p == 2
        return LinearAlgebra.norm2(x)
    elseif p == Inf
        return maximum(abs(x))
    elseif p > 1
        # TODO: allow tolerance in the rationalize step
        return rationalnorm(x, rationalize(Int, float(p)))
    else
        error("vector p-norms not defined for p < 1")
    end
end
