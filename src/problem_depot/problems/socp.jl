# Copyright (c) 2014: Madeleine Udell and contributors
#
# Use of this source code is governed by a BSD-style license that can be found
# in the LICENSE file or at https://opensource.org/license/bsd-2-clause

@add_problem socp function socp_dual_norm_2_atom(
    handle_problem!,
    ::Val{test},
    atol,
    rtol,
    ::Type{T},
) where {T,test}
    x = Variable(2, 1)
    A = [1 2; 2 1; 3 4]
    b = [2; 3; 4]
    p = minimize(norm2(A * x + b); numeric_type = T)

    if test
        @test problem_vexity(p) == ConvexVexity()
    end
    handle_problem!(p)
    if test
        @test p.optval ≈ 0.64888 atol = atol rtol = rtol
        @test evaluate(norm2(A * x + b)) ≈ 0.64888 atol = atol rtol = rtol
    end

    x = Variable(2, 1)
    A = [1 2; 2 1; 3 4]
    b = [2; 3; 4]
    lambda = 1
    p = minimize(norm2(A * x + b) + lambda * norm2(x), x >= 1; numeric_type = T)

    if test
        @test problem_vexity(p) == ConvexVexity()
    end
    handle_problem!(p)
    if test
        @test p.optval ≈ 14.9049 atol = atol rtol = rtol
        @test evaluate(norm2(A * x + b) + lambda * norm2(x)) ≈ 14.9049 atol =
            atol rtol = rtol
        @test p.constraints[1].dual ≈ [4.4134, 5.1546] atol = atol rtol = rtol
    end

    x = Variable(2)

    p = minimize(
        norm2([x[1] + 2x[2] + 2; 2x[1] + x[2] + 3; 3x[1] + 4x[2] + 4]) +
        lambda * norm2(x),
        x >= 1;
        numeric_type = T,
    )

    if test
        @test problem_vexity(p) == ConvexVexity()
    end

    handle_problem!(p)
    if test
        @test p.optval ≈ 14.9049 atol = atol rtol = rtol
        @test evaluate(norm2(A * x + b) + lambda * norm2(x)) ≈ 14.9049 atol =
            atol rtol = rtol
        @test p.constraints[1].dual ≈ [4.4134, 5.1546] atol = atol rtol = rtol
    end

    x = Variable(2, 1)
    A = [1 2; 2 1; 3 4]
    b = [2; 3; 4]
    lambda = 1
    p = minimize(
        norm2(A * x + b) + lambda * norm(x, 1),
        x >= 1;
        numeric_type = T,
    )

    if test
        @test problem_vexity(p) == ConvexVexity()
    end
    handle_problem!(p)
    if test
        @test p.optval ≈ 15.4907 atol = atol rtol = rtol
        @test evaluate(norm2(A * x + b) + lambda * norm(x, 1)) ≈ 15.4907 atol =
            atol rtol = rtol
        @test p.constraints[1].dual ≈ [4.7062, 5.4475] atol = atol rtol = rtol
    end
end

@add_problem socp function socp_dual_frobenius_norm_atom(
    handle_problem!,
    ::Val{test},
    atol,
    rtol,
    ::Type{T},
) where {T,test}
    m = Variable(4, 5)
    c = [m[3, 3] == 4, m >= 1]
    p = minimize(norm(m, 2), c; numeric_type = T)

    if test
        @test problem_vexity(p) == ConvexVexity()
    end
    handle_problem!(p)
    if test
        @test p.optval ≈ sqrt(35) atol = atol rtol = rtol
        @test evaluate(norm(m, 2)) ≈ sqrt(35) atol = atol rtol = rtol
        @test p.constraints[1].dual ≈ 0.6761 atol = atol rtol = rtol
        dual = 0.1690 .* ones(4, 5)
        dual[3, 3] = 0
        @test p.constraints[2].dual ≈ dual atol = atol rtol = rtol
    end
end

@add_problem socp function socp_quad_over_lin_atom(
    handle_problem!,
    ::Val{test},
    atol,
    rtol,
    ::Type{T},
) where {T,test}
    x = Variable(3)
    A = [2 -3 5; -2 9 -3; 5 -8 3]
    b = [-3; 9; 5]
    c = [3, 2, 4]
    d = -3
    p = minimize(quadoverlin(A * x + b, dot(c, x) + d); numeric_type = T)

    if test
        @test problem_vexity(p) == ConvexVexity()
    end
    handle_problem!(p)
    if test
        @test p.optval ≈ 17.7831 atol = atol rtol = rtol
        @test evaluate(quadoverlin(A * x + b, dot(c, x) + d)) ≈ 17.7831 atol =
            atol rtol = rtol
    end
end

@add_problem socp function socp_sum_squares_atom(
    handle_problem!,
    ::Val{test},
    atol,
    rtol,
    ::Type{T},
) where {T,test}
    x = Variable(2, 1)
    A = [1 2; 2 1; 3 4]
    b = [2; 3; 4]
    p = minimize(sumsquares(A * x + b); numeric_type = T)

    if test
        @test problem_vexity(p) == ConvexVexity()
    end
    handle_problem!(p)
    if test
        @test p.optval ≈ 0.42105 atol = atol rtol = rtol
        @test (evaluate(sumsquares(A * x + b)))[1] ≈ 0.42105 atol = atol rtol =
            rtol
    end
end

@add_problem socp function socp_square_atom(
    handle_problem!,
    ::Val{test},
    atol,
    rtol,
    ::Type{T},
) where {T,test}
    x = Variable(2, 1)
    A = [1 2; 2 1; 3 4]
    b = [2; 3; 4]
    p = minimize(sum(square(A * x + b)); numeric_type = T)

    if test
        @test problem_vexity(p) == ConvexVexity()
    end
    handle_problem!(p)
    if test
        @test p.optval ≈ 0.42105 atol = atol rtol = rtol
        @test evaluate(sum(square(A * x + b))) ≈ 0.42105 atol = atol rtol = rtol
    end

    x = Variable(2, 1)
    A = [1 2; 2 1; 3 4]
    b = [2; 3; 4]
    expr = A * x + b
    # `literal_pow` case:
    p = minimize(sum(expr .^ 2); numeric_type = T) # elementwise ^
    if test
        @test problem_vexity(p) == ConvexVexity()
    end
    handle_problem!(p)
    if test
        @test p.optval ≈ 0.42105 atol = atol rtol = rtol
        @test evaluate(sum(broadcast(^, expr, 2))) ≈ 0.42105 atol = atol rtol =
            rtol
    end

    # Test non-literal case:
    k = 2
    p = minimize(sum(expr .^ k); numeric_type = T) # elementwise ^
    if test
        @test problem_vexity(p) == ConvexVexity()
    end
    handle_problem!(p)
    if test
        @test p.optval ≈ 0.42105 atol = atol rtol = rtol
        @test evaluate(sum(broadcast(^, expr, 2))) ≈ 0.42105 atol = atol rtol =
            rtol
    end

    p = minimize(sum(expr .* expr); numeric_type = T) # elementwise *
    if test
        @test problem_vexity(p) == ConvexVexity()
    end
    handle_problem!(p)
    if test
        @test p.optval ≈ 0.42105 atol = atol rtol = rtol
        @test evaluate(sum(expr .* expr)) ≈ 0.42105 atol = atol rtol = rtol
    end
end

@add_problem socp function socp_inv_pos_atom(
    handle_problem!,
    ::Val{test},
    atol,
    rtol,
    ::Type{T},
) where {T,test}
    x = Variable(4)
    p = minimize(
        sum(invpos(x)),
        invpos(x) <= 2,
        x >= 1,
        x == 2,
        2 == x;
        numeric_type = T,
    )

    if test
        @test problem_vexity(p) == ConvexVexity()
    end
    handle_problem!(p)
    if test
        @test p.optval ≈ 2 atol = atol rtol = rtol
        @test evaluate(sum(invpos(x))) ≈ 2 atol = atol rtol = rtol
    end

    x = Variable(3)
    p = minimize(sum([3, 6, 9] ./ x), x <= 3; numeric_type = T)

    handle_problem!(p)
    if test
        @test evaluate(x) ≈ fill(3.0, (3, 1)) atol = atol rtol = rtol
        @test p.optval ≈ 6 atol = atol rtol = rtol
        @test evaluate(sum([3, 6, 9] ./ x)) ≈ 6 atol = atol rtol = rtol
    end

    x = Variable()
    p = minimize(sum([3, 6, 9] / x), x <= 3; numeric_type = T)

    handle_problem!(p)
    if test
        @test evaluate(x) ≈ 3 atol = atol rtol = rtol
        @test p.optval ≈ 6 atol = atol rtol = rtol
        @test evaluate(sum([3, 6, 9] / x)) ≈ 6 atol = atol rtol = rtol
    end
end

@add_problem socp function socp_geo_mean_atom(
    handle_problem!,
    ::Val{test},
    atol,
    rtol,
    ::Type{T},
) where {T,test}
    x = Variable(2)
    y = Variable(2)
    p = minimize(geomean(x, y), x >= 1, y >= 2; numeric_type = T)

    # not DCP compliant
    if test
        @test problem_vexity(p) == ConcaveVexity()
    end
    p = maximize(geomean(x, y), 1 <= x, x <= 2, y <= 2; numeric_type = T)

    # Just gave it a vector as an objective, not okay
    if test
        @test_throws Exception handle_problem!(p)
    end

    p = maximize(sum(geomean(x, y)), 1 <= x, x <= 2, y <= 2; numeric_type = T)

    handle_problem!(p)
    if test
        @test p.optval ≈ 4 atol = atol rtol = rtol
        @test evaluate(sum(geomean(x, y))) ≈ 4 atol = atol rtol = rtol
    end

    # 3 arg
    z = Variable(2)
    p = maximize(
        sum(geomean(x, y, z)),
        1 <= x,
        x <= 2,
        y <= 2,
        z <= 2;
        numeric_type = T,
    )
    handle_problem!(p)
    if test
        @test p.optval ≈ 4 atol = atol rtol = rtol
        @test evaluate(sum(geomean(x, y, z))) ≈ 4 atol = atol rtol = rtol
    end

    p = maximize(
        sum(geomean(x, y, 4 * ones(2))),
        1 <= x,
        x <= 2,
        y <= 2;
        numeric_type = T,
    )
    handle_problem!(p)
    if test
        @test p.optval ≈ 2 * 4^(2 / 3) atol = atol rtol = rtol
        @test evaluate(sum(geomean(x, y, 4 * ones(2)))) ≈ 2 * 4^(2 / 3) atol =
            atol rtol = rtol
    end
end

@add_problem socp function socp_sqrt_atom(
    handle_problem!,
    ::Val{test},
    atol,
    rtol,
    ::Type{T},
) where {T,test}
    x = Variable()
    return p = maximize(sqrt(x), 1 >= x; numeric_type = T)
end

@add_problem socp function socp_quad_form_atom(
    handle_problem!,
    ::Val{test},
    atol,
    rtol,
    ::Type{T},
) where {T,test}
    x = Variable(3, 1)
    A = [0.8608 0.3131 0.5458; 0.3131 0.8584 0.5836; 0.5458 0.5836 1.5422]
    p = minimize(quadform(x, A), [x >= 1]; numeric_type = T)

    if test
        @test problem_vexity(p) == ConvexVexity()
    end
    handle_problem!(p)
    if test
        @test p.optval ≈ 6.1464 atol = atol rtol = rtol
        @test evaluate(quadform(x, A)) ≈ 6.1464 atol = atol rtol = rtol
    end

    x = Variable(3, 1)
    A =
        -1.0 *
        [0.8608 0.3131 0.5458; 0.3131 0.8584 0.5836; 0.5458 0.5836 1.5422]
    c = [3 2 4]
    p = maximize(c * x, [quadform(x, A) >= -1]; numeric_type = T)

    if test
        @test problem_vexity(p) == ConvexVexity()
    end
    handle_problem!(p)
    if test
        @test p.optval ≈ 3.7713 atol = atol rtol = rtol
        @test evaluate(quadform(x, A)) ≈ -1 atol = atol rtol = rtol
    end

    # https://github.com/jump-dev/Convex.jl/issues/398
    n = 3
    b = randn(n)
    x = Variable(n)
    H = Semidefinite(n)
    Hval = randn(n, n)
    Hval .= Hval' * Hval + 10 * diagm(0 => ones(n)) # symmetric positive definite
    fix!(H, Hval)
    p = minimize(x'b + quadform(x, evaluate(H)), [x >= 0]; numeric_type = T)
    handle_problem!(p)

    p2 = minimize(x'b + quadform(x, Hval), [x >= 0]; numeric_type = T)
    handle_problem!(p2)

    if test
        @test p.optval ≈ p2.optval atol = atol rtol = rtol
        @test evaluate(H) ≈ Hval atol = atol rtol = rtol
    end

    # https://github.com/jump-dev/Convex.jl/pull/444
    x = Variable(3)
    a, b, c, d, e, f = rand(6)
    M = [
        3 a-b*im c-d*im
        a+b*im 3 e-f*im
        c+d*im e+f*im 3
    ]
    y = rand(3)
    p = minimize(quadform(x - y, M); numeric_type = T)
    handle_problem!(p)
    if test
        @test evaluate(x) ≈ y atol = atol rtol = rtol
    end

    # https://github.com/jump-dev/Convex.jl/issues/446
    x = Variable(3)
    A = zeros(3, 3)
    A[1, 1] = 1.0
    p = minimize(
        quadform(x, A) - 2 * x[2] - x[3],
        norm(x, 1) <= 1;
        numeric_type = T,
    )
    handle_problem!(p)
    if test
        @test evaluate(x) ≈ [0.0, 1.0, 0.0] atol = atol rtol = rtol
    end
end

@add_problem socp function socp_huber_atom(
    handle_problem!,
    ::Val{test},
    atol,
    rtol,
    ::Type{T},
) where {T,test}
    x = Variable(3)
    p = minimize(sum(huber(x, 1)), x >= 2; numeric_type = T)

    if test
        @test problem_vexity(p) == ConvexVexity()
    end
    handle_problem!(p)
    if test
        @test p.optval ≈ 9 atol = atol rtol = rtol
        @test evaluate(sum(huber(x, 1))) ≈ 9 atol = atol rtol = rtol
    end
end

@add_problem socp function socp_rational_norm_atom(
    handle_problem!,
    ::Val{test},
    atol,
    rtol,
    ::Type{T},
) where {T,test}
    A = [1 2 3; -1 2 3]
    b = A * ones(3)
    x = Variable(3)
    p = minimize(norm(x, 4.5), [A * x == b]; numeric_type = T)

    if test
        @test problem_vexity(p) == ConvexVexity()
    end
    # Solution is approximately x = [1, .93138, 1.04575]
    handle_problem!(p)
    if test
        @test p.optval ≈ 1.2717 atol = atol rtol = rtol
        @test evaluate(norm(x, 4.5)) ≈ 1.2717 atol = atol rtol = rtol
    end
end

@add_problem socp function socp_rational_norm_dual_norm(
    handle_problem!,
    ::Val{test},
    atol,
    rtol,
    ::Type{T},
) where {T,test}
    v = [0.463339, 0.0216084, -2.07914, 0.99581, 0.889391]
    x = Variable(5)
    q = 1.379  # q norm constraint that generates many inequalities
    qs = q / (q - 1)  # Conjugate to q
    p = minimize(x' * v; numeric_type = T)

    p.constraints += (norm(x, q) <= 1)
    if test
        @test problem_vexity(p) == ConvexVexity()
    end
    handle_problem!(p) # Solution is -norm(v, q / (q - 1))
    if test
        @test p.optval ≈ -2.144087 atol = atol rtol = rtol
        @test sum(evaluate(x' * v)) ≈ -2.144087 atol = atol rtol = rtol
        @test evaluate(norm(x, q)) ≈ 1 atol = atol rtol = rtol
        @test sum(evaluate(x' * v)) ≈ -(sum(abs.(v) .^ qs)^(1 / qs)) atol = atol rtol =
            rtol
    end
end

@add_problem socp function socp_rational_norm_atom_sum(
    handle_problem!,
    ::Val{test},
    atol,
    rtol,
    ::Type{T},
) where {T,test}
    A = [
        -0.719255 -0.229089
        -1.33632 -1.37121
        0.703447 -1.4482
    ]
    b = [-1.82041, -1.67516, -0.866884]
    q = 1.5
    xvar = Variable(2)
    p = minimize(
        0.5 * sumsquares(xvar) + norm(A * xvar - b, q);
        numeric_type = T,
    )

    if test
        @test problem_vexity(p) == ConvexVexity()
    end
    handle_problem!(p)

    if test
        # Compute gradient, check it is zero(ish)
        x_opt = evaluate(xvar)
        margins = A * x_opt - b
        qs = q / (q - 1)  # Conjugate
        denom = sum(abs.(margins) .^ q)^(1 / qs)
        g = x_opt + A' * (abs.(margins) .^ (q - 1) .* sign.(margins)) / denom
        @test p.optval ≈ 1.7227 atol = atol rtol = rtol
        @test norm(g, 2)^2 ≈ 0 atol = atol rtol = rtol
    end
end

@add_problem socp function socp_norm_consistent_with_Base_for_matrix_variables(
    handle_problem!,
    ::Val{test},
    atol,
    rtol,
    ::Type{T},
) where {T,test}
    A = randn(4, 4)
    x = Variable(4, 4)
    set_value!(x, A)
    # Matrix norm
    if test
        @test evaluate(LinearAlgebra.opnorm(x)) ≈ LinearAlgebra.opnorm(A) atol =
            atol rtol = rtol
        @test evaluate(LinearAlgebra.opnorm(x, 1)) ≈ LinearAlgebra.opnorm(A, 1) atol =
            atol rtol = rtol
        @test evaluate(LinearAlgebra.opnorm(x, 2)) ≈ LinearAlgebra.opnorm(A, 2) atol =
            atol rtol = rtol
        @test evaluate(LinearAlgebra.opnorm(x, Inf)) ≈
              LinearAlgebra.opnorm(A, Inf) atol = atol rtol = rtol
    end
    # Vector norm
    if test
        @test evaluate(norm(x, 1)) ≈ norm(vec(A), 1) atol = atol rtol = rtol
        @test evaluate(norm(x, 2)) ≈ norm(vec(A), 2) atol = atol rtol = rtol
        @test evaluate(norm(x, 7)) ≈ norm(vec(A), 7) atol = atol rtol = rtol
        @test evaluate(norm(x, Inf)) ≈ norm(vec(A), Inf) atol = atol rtol = rtol
    end
end

@add_problem socp function socp_fix_and_free_addition(
    handle_problem!,
    ::Val{test},
    atol,
    rtol,
    ::Type{T},
) where {T,test}
    x = Variable()
    y = Variable()

    p = minimize(x + y, x >= 0, y >= 0; numeric_type = T)

    handle_problem!(p)
    if test
        @test p.optval ≈ 0 atol = atol rtol = rtol
    end

    set_value!(y, 4)
    fix!(y)
    handle_problem!(p)
    if test
        @test p.optval ≈ 4 atol = atol rtol = rtol
    end

    free!(y)
    handle_problem!(p)
    if test
        @test p.optval ≈ 0 atol = atol rtol = rtol
    end
end

@add_problem socp function socp_fix_multiplication(
    handle_problem!,
    ::Val{test},
    atol,
    rtol,
    ::Type{T},
) where {T,test}
    a = [1, 2, 3, 2, 1]
    x = Variable(length(a))
    gamma = Variable(Positive())
    fix!(gamma, 0.7)

    p = minimize(
        norm(x - a) + gamma * norm(x[1:end-1] - x[2:end]);
        numeric_type = T,
    )

    handle_problem!(p)
    if test
        o1 = p.optval
        # x should be very close to a
        @test o1 ≈ 0.7 * norm(a[1:end-1] - a[2:end]) atol = atol rtol = rtol
    end
    # increase regularization
    fix!(gamma, 1.0)
    handle_problem!(p)

    if test
        o2 = p.optval
        # x should be very close to mean(a)
        @test o2 ≈ norm(a .- mean(a)) atol = atol rtol = rtol
    end

    if test
        @test o1 <= o2
    end
end

@add_problem socp function socp_dual_minimal_norm_solutions(
    handle_problem!,
    ::Val{test},
    atol,
    rtol,
    ::Type{T},
) where {T,test}
    x = Variable(2)
    A = [1 2; 2 4]
    b = [3, 6]
    p = minimize(norm(x, 1), A * x == b; numeric_type = T)

    if test
        @test problem_vexity(p) == ConvexVexity()
    end

    handle_problem!(p)
    if test
        @test p.optval ≈ 1.5 atol = atol rtol = rtol
        @test evaluate(x) ≈ [0, 1.5] atol = atol rtol = rtol
        @test evaluate(norm(x, 1)) ≈ p.optval atol = atol rtol = rtol
        @test dot(b, p.constraints[1].dual) ≈ p.optval atol = atol rtol = rtol
    end

    x = Variable(2)
    A = [1 2; 2 4]
    b = [3, 6]
    p = minimize(norm(x, 2), A * x == b; numeric_type = T)

    test && @test problem_vexity(p) == ConvexVexity()

    handle_problem!(p)

    if test
        @test p.optval ≈ 3 / sqrt(5) atol = atol rtol = rtol
        @test evaluate(x) ≈ [3 / 5, 6 / 5] atol = atol rtol = rtol
        @test evaluate(norm(x, 2)) ≈ p.optval atol = atol rtol = rtol
        @test dot(b, p.constraints[1].dual) ≈ p.optval atol = atol rtol = rtol
    end

    x = Variable(2)
    A = [1 2; 2 4]
    b = [3, 6]
    p = minimize(norm(x, Inf), A * x == b; numeric_type = T)

    test && @test problem_vexity(p) == ConvexVexity()

    handle_problem!(p)

    if test
        @test p.optval ≈ 1.0 atol = atol rtol = rtol
        @test evaluate(x) ≈ [1, 1] atol = atol rtol = rtol
        @test evaluate(norm(x, Inf)) ≈ p.optval atol = atol rtol = rtol
        @test dot(b, p.constraints[1].dual) ≈ p.optval atol = atol rtol = rtol
    end
end
