# Copyright (c) 2014: Madeleine Udell and contributors
#
# Use of this source code is governed by a BSD-style license that can be found
# in the LICENSE file or at https://opensource.org/license/bsd-2-clause

LinearAlgebra.tr(e::AbstractExpr) = sum(LinearAlgebra.diag(e))
