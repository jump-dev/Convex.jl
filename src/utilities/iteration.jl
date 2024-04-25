# Copyright (c) 2014: Madeleine Udell and contributors
#
# Use of this source code is governed by a BSD-style license that can be found
# in the LICENSE file or at https://opensource.org/license/bsd-2-clause

function Base.iterate(x::AbstractVariable, s = 0)
    return s >= length(x) ? nothing : (x[s+1], s + 1)
end
