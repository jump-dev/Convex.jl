# Copyright (c) 2014: Madeleine Udell and contributors
#
# Use of this source code is governed by a BSD-style license that can be found
# in the LICENSE file or at https://opensource.org/license/bsd-2-clause

# Using _vec avoids broadcast-specific behaviors that we want to avoid,
# such as extending singleton dimensions. We need to ensure that the inputs have
# the same length, which broadcast will check for us if both inputs are vectors.
_vec(x::AbstractExpr) = size(x, 1) > 1 && size(x, 2) > 1 ? vec(x) : x

_vec(x::AbstractMatrix) = vec(x)

_vec(x::Any) = x

_expr_vec(x) = convert(AbstractExpr, _vec(x))

_dot(x, y) = sum(broadcast(*, conj(_expr_vec(x)), _expr_vec(y)))

"""
    LinearAlgebra.dot(x::Convex.AbstractExpr, y::Convex.AbstractExpr)

The dot product \$x \\cdot y\$. If `x` is complex, it is conjugated.

## Examples

```jldoctest; filter=r"id: [0-9]+…[0-9]+"
julia> x = ComplexVariable(2);

julia> y = [1, 2];

julia> atom = dot(x, y)
sum (affine; complex)
└─ .* (affine; complex)
   ├─ conj (affine; complex)
   │  └─ 2-element complex variable (id: 133…443)
   └─ [1; 2;;]

julia> size(atom)
(1, 1)
```
"""
LinearAlgebra.dot(x::AbstractExpr, y::AbstractExpr) = _dot(x, y)

LinearAlgebra.dot(x::Value, y::AbstractExpr) = _dot(x, y)

LinearAlgebra.dot(x::AbstractExpr, y::Value) = _dot(x, y)
