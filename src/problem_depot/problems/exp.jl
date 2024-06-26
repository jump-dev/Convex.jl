# Copyright (c) 2014: Madeleine Udell and contributors
#
# Use of this source code is governed by a BSD-style license that can be found
# in the LICENSE file or at https://opensource.org/license/bsd-2-clause

@add_problem exp function exp_exp_atom(
    handle_problem!,
    ::Val{test},
    atol,
    rtol,
    ::Type{T},
) where {T,test}
    y = Variable()
    p = minimize(exp(y), y >= 0; numeric_type = T)

    if test
        @test problem_vexity(p) == ConvexVexity()
    end
    handle_problem!(p)
    if test
        @test p.optval ≈ 1 atol = atol rtol = rtol
        @test evaluate(exp(y)) ≈ 1 atol = atol rtol = rtol
    end

    # Test for constant `exp` (#613)
    y = Variable()
    x = constant([1, 2, 3])
    p = minimize(sum(exp(x)) + y, y >= 0; numeric_type = T)
    if test
        @test problem_vexity(p) == ConvexVexity()
    end
    handle_problem!(p)
    if test
        @test p.optval ≈ sum(exp.([1, 2, 3])) atol = atol rtol = rtol
        @test evaluate(y) ≈ 0 atol = atol
    end

    y = Variable()
    p = minimize(exp(y), y >= 1; numeric_type = T)

    if test
        @test problem_vexity(p) == ConvexVexity()
    end
    handle_problem!(p)
    if test
        @test p.optval ≈ exp(1) atol = atol rtol = rtol
        @test evaluate(exp(y)) ≈ exp(1) atol = atol rtol = rtol
    end

    y = Variable(5)
    p = minimize(sum(exp(y)), y >= 0; numeric_type = T)

    if test
        @test problem_vexity(p) == ConvexVexity()
    end
    handle_problem!(p)
    if test
        @test p.optval ≈ 5 atol = atol rtol = rtol
        @test evaluate(sum(exp(y))) ≈ 5 atol = atol rtol = rtol
    end

    y = Variable(5)
    p = minimize(sum(exp(y)), y >= 0; numeric_type = T)

    if test
        @test problem_vexity(p) == ConvexVexity()
    end
    handle_problem!(p)
    if test
        @test p.optval ≈ 5 atol = atol rtol = rtol
    end
end

@add_problem exp function exp_log_atom(
    handle_problem!,
    ::Val{test},
    atol,
    rtol,
    ::Type{T},
) where {T,test}
    y = Variable()
    p = maximize(log(y), y <= 1; numeric_type = T)

    if test
        @test problem_vexity(p) == ConvexVexity()
    end
    handle_problem!(p)
    if test
        @test p.optval ≈ 0 atol = atol rtol = rtol
    end

    y = Variable()
    p = maximize(log(y), y <= 2; numeric_type = T)

    if test
        @test problem_vexity(p) == ConvexVexity()
    end
    handle_problem!(p)
    if test
        @test p.optval ≈ log(2) atol = atol rtol = rtol
    end

    y = Variable()
    p = maximize(log(y), [y <= 2, exp(y) <= 10]; numeric_type = T)

    if test
        @test problem_vexity(p) == ConvexVexity()
    end
    handle_problem!(p)
    if test
        @test p.optval ≈ log(2) atol = atol rtol = rtol
    end
end

@add_problem exp function exp_log_sum_exp_atom(
    handle_problem!,
    ::Val{test},
    atol,
    rtol,
    ::Type{T},
) where {T,test}
    y = Variable(5)
    p = minimize(logsumexp(y), y >= 1; numeric_type = T)
    if test
        @test problem_vexity(p) == ConvexVexity()
    end
    handle_problem!(p)
    if test
        @test p.optval ≈ log(exp(1) * 5) atol = atol rtol = rtol
    end
    y = Variable(5, 2)
    p = minimize(
        sum(Convex.logsumexp(y; dims = 1)),
        y[:, 1] >= 1,
        y[:, 2] >= 2;
        numeric_type = T,
    )
    handle_problem!(p)
    if test
        @test evaluate(y[:, 1]) ≈ ones(5) atol = atol rtol = rtol
        @test evaluate(y[:, 2]) ≈ 2 * ones(5) atol = atol rtol = rtol
        @test ≈(
            p.optval,
            log(exp(1) * 5) + log(exp(2) * 5);
            atol = atol,
            rtol = rtol,
        )
    end
    p = minimize(logsumexp(y), y[:, 1] >= 1, y[:, 2] >= 2; numeric_type = T)
    handle_problem!(p)
    if test
        @test evaluate(y[:, 1]) ≈ ones(5) atol = atol rtol = rtol
        @test evaluate(y[:, 2]) ≈ 2 * ones(5) atol = atol rtol = rtol
        @test p.optval ≈ log(exp(1) * 5 + exp(2) * 5) atol = atol rtol = rtol
    end

    x = Variable(2, 3)
    v = Convex.logsumexp(x; dims = 1)
    p = minimize(sum(v), x >= [1 2 3; 4 5 6]; numeric_type = T)
    handle_problem!(p)
    if test
        @test evaluate(x) ≈ [1 2 3; 4 5 6] atol = atol rtol = rtol
        @test ≈(
            evaluate(v),
            log.(sum(exp, evaluate(x); dims = 1));
            atol = atol,
            rtol = rtol,
        )
        @test vexity(v) == Convex.ConvexVexity()
    end
    return
end

@add_problem exp function exp_logistic_loss_atom(
    handle_problem!,
    ::Val{test},
    atol,
    rtol,
    ::Type{T},
) where {T,test}
    y = Variable(5)
    p = minimize(logisticloss(y), y >= 1; numeric_type = T)

    if test
        @test problem_vexity(p) == ConvexVexity()
    end
    handle_problem!(p)
    if test
        @test p.optval ≈ log(exp(1) + 1) * 5 atol = atol rtol = rtol
    end
end

@add_problem exp function exp_entropy_atom(
    handle_problem!,
    ::Val{test},
    atol,
    rtol,
    ::Type{T},
) where {T,test}
    y = Variable(5, Positive())
    p = maximize(entropy(y), sum(y) <= 1; numeric_type = T)

    if test
        @test problem_vexity(p) == ConvexVexity()
    end
    handle_problem!(p)
    if test
        @test p.optval ≈ -(log(1 / 5)) atol = atol rtol = rtol
    end
end

@add_problem exp function exp_relative_entropy_atom(
    handle_problem!,
    ::Val{test},
    atol,
    rtol,
    ::Type{T},
) where {T,test}
    x = Variable(1)
    y = Variable(1)
    # x log (x/y)
    p = minimize(relative_entropy(x, y), y == 1, x >= 2; numeric_type = T)

    if test
        @test problem_vexity(p) == ConvexVexity()
    end
    handle_problem!(p)
    if test
        @test p.optval ≈ 2 * log(2) atol = atol rtol = rtol
    end
end

@add_problem exp function exp_log_perspective_atom(
    handle_problem!,
    ::Val{test},
    atol,
    rtol,
    ::Type{T},
) where {T,test}
    x = Variable(1)
    y = Variable(1)
    # y log (x/y)
    p = maximize(log_perspective(x, y), y == 5, x <= 10; numeric_type = T)

    if test
        @test problem_vexity(p) == ConvexVexity()
    end
    handle_problem!(p)
    if test
        @test p.optval ≈ 5 * log(2) atol = atol rtol = rtol
    end
end
