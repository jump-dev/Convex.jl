# Copyright (c) 2014: Madeleine Udell and contributors
#
# Use of this source code is governed by a BSD-style license that can be found
# in the LICENSE file or at https://opensource.org/license/bsd-2-clause

function inner_product(x::AbstractExpr, y::AbstractExpr)
    if !(x.size == y.size && x.size[1] == x.size[2])
        error("arguments must be square matrices of the same dimension")
    end
    return real(LinearAlgebra.tr(x' * y))
end

inner_product(x::Value, y::AbstractExpr) = inner_product(constant(x), y)

inner_product(x::AbstractExpr, y::Value) = inner_product(x, constant(y))
