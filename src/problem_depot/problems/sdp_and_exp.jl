# Copyright (c) 2014: Madeleine Udell and contributors
#
# Use of this source code is governed by a BSD-style license that can be found
# in the LICENSE file or at https://opensource.org/license/bsd-2-clause

@add_problem sdp_and_exp function sdp_and_exp_log_det_atom(
    handle_problem!,
    ::Val{test},
    atol,
    rtol,
    ::Type{T},
) where {T,test}
    x = Variable(2, 2)
    p = maximize(logdet(x), [x[1, 1] == 1, x[2, 2] == 1]; numeric_type = T)

    if test
        @test problem_vexity(p) == ConvexVexity()
    end
    handle_problem!(p)
    if test
        @test p.optval ≈ 0 atol = 10atol rtol = rtol
        @test evaluate(logdet(x)) ≈ 0 atol = 10atol rtol = rtol
    end
end
