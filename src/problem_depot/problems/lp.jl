@add_problem lp function lp_dual_abs_atom(
    handle_problem!,
    ::Val{test},
    atol,
    rtol,
    ::Type{T},
) where {T,test}
    x = Variable()
    p = minimize(abs(x), x <= -1; numeric_type = T)

    if test
        @test problem_vexity(p) == ConvexVexity()
    end
    handle_problem!(p)
    if test
        @test objective_value(p) ≈ 1 atol = atol rtol = rtol
        @test evaluate(abs(x)) ≈ 1 atol = atol rtol = rtol
        @test p.constraints[1].dual ≈ -1 atol = atol rtol = rtol
    end

    x = Variable(2, 2)
    p = minimize(
        sum(abs(x)),
        x[2, 2] >= 1,
        x[1, 1] >= 1,
        x >= 0;
        numeric_type = T,
    )

    if test
        @test problem_vexity(p) == ConvexVexity()
    end
    handle_problem!(p)
    if test
        @test objective_value(p) ≈ 2 atol = atol rtol = rtol
        @test evaluate(sum(abs(x))) ≈ 2 atol = atol rtol = rtol
        @test p.constraints[1].dual ≈ 1 atol = atol rtol = rtol
        @test p.constraints[2].dual ≈ 1 atol = atol rtol = rtol
        @test p.constraints[3].dual[1, 1] ≈ 0 atol = atol rtol = rtol
        @test p.constraints[3].dual[2, 2] ≈ 0 atol = atol rtol = rtol
        @test p.constraints[3].dual[1, 2] ≈ p.constraints[3].dual[2, 1] atol =
            atol rtol = rtol
        @test p.constraints[3].dual[1, 2] <= 1 + atol
        @test p.constraints[3].dual[2, 1] <= 1 + atol
        @test p.constraints[3].dual[1, 2] >= -atol
        @test p.constraints[3].dual[2, 1] >= -atol
    end
end

@add_problem lp function lp_maximum_atom(
    handle_problem!,
    ::Val{test},
    atol,
    rtol,
    ::Type{T},
) where {T,test}
    x = Variable(10)
    a = shuffle(collect(0.1:0.1:1.0))
    p = minimize(maximum(x), x >= a; numeric_type = T)

    if test
        @test problem_vexity(p) == ConvexVexity()
    end
    handle_problem!(p)
    if test
        @test objective_value(p) ≈ maximum(a) atol = atol rtol = rtol
        @test evaluate(maximum(x)) ≈ maximum(a) atol = atol rtol = rtol
    end
end

@add_problem lp function lp_minimum_atom(
    handle_problem!,
    ::Val{test},
    atol,
    rtol,
    ::Type{T},
) where {T,test}
    x = Variable(1)
    a = reshape(shuffle(collect(0.01:0.01:1.0)), (10, 10))
    p = maximize(minimum(x), x <= a; numeric_type = T)

    if test
        @test problem_vexity(p) == ConvexVexity()
    end
    handle_problem!(p)
    if test
        @test objective_value(p) ≈ minimum(a) atol = atol rtol = rtol
        @test evaluate(minimum(x)) ≈ minimum(a) atol = atol rtol = rtol
    end

    x = Variable(4, 4)
    y = Variable(4, 6)
    z = Variable(1)
    c = ones(4, 1)
    d = fill(2.0, (6, 1))
    constraints = [[x y] <= 2, z <= 0, z <= x, 2z >= -1]
    objective = sum(x + z * ones(4, 4)) + minimum(y) + c' * y * d
    p = maximize(objective, constraints; numeric_type = T)

    if test
        @test problem_vexity(p) == ConvexVexity()
    end
    handle_problem!(p)
    if test
        @test objective_value(p) ≈ 130 atol = atol rtol = rtol
        @test (evaluate(objective))[1] ≈ 130 atol = atol rtol = rtol
    end
end

@add_problem lp function lp_max_atom(
    handle_problem!,
    ::Val{test},
    atol,
    rtol,
    ::Type{T},
) where {T,test}
    x = Variable(10, 10)
    y = Variable(10, 10)
    a = reshape(shuffle(collect(0.01:0.01:1.0)), (10, 10))
    b = reshape(shuffle(collect(0.01:0.01:1.0)), (10, 10))
    p = minimize(maximum(max(x, y)), [x >= a, y >= b]; numeric_type = T)

    if test
        @test problem_vexity(p) == ConvexVexity()
    end
    handle_problem!(p)
    max_a = maximum(a)
    max_b = maximum(b)
    if test
        @test objective_value(p) ≈ max(max_a, max_b) atol = 10atol atol = atol rtol =
            rtol
        @test evaluate(maximum(max(x, y))) ≈ max(max_a, max_b) atol = 10atol atol =
            atol rtol = rtol
    end
end

@add_problem lp function lp_min_atom(
    handle_problem!,
    ::Val{test},
    atol,
    rtol,
    ::Type{T},
) where {T,test}
    x = Variable(10, 10)
    y = Variable(10, 10)
    a = reshape(shuffle(collect(0.01:0.01:1.0)), (10, 10))
    b = reshape(shuffle(collect(0.01:0.01:1.0)), (10, 10))
    p = maximize(minimum(min(x, y)), [x <= a, y <= b]; numeric_type = T)

    if test
        @test problem_vexity(p) == ConvexVexity()
    end
    handle_problem!(p)
    min_a = minimum(a)
    min_b = minimum(b)
    if test
        @test objective_value(p) ≈ min(min_a, min_b) atol = 10atol atol = atol rtol =
            rtol
        @test evaluate(minimum(min(x, y))) ≈ min(min_a, min_b) atol = 10atol atol =
            atol rtol = rtol
    end
end

@add_problem lp function lp_pos_atom(
    handle_problem!,
    ::Val{test},
    atol,
    rtol,
    ::Type{T},
) where {T,test}
    x = Variable(3)
    a = [-2; 1; 2]
    p = minimize(sum(pos(x)), [x >= a, x <= 2]; numeric_type = T)

    if test
        @test problem_vexity(p) == ConvexVexity()
    end
    handle_problem!(p)
    if test
        @test objective_value(p) ≈ 3 atol = atol rtol = rtol
        @test evaluate(sum(pos(x))) ≈ 3 atol = atol rtol = rtol
    end
end

@add_problem lp function lp_neg_atom(
    handle_problem!,
    ::Val{test},
    atol,
    rtol,
    ::Type{T},
) where {T,test}
    x = Variable(3)
    p = minimize(1, [x >= -2, x <= -2, neg(x) <= 3]; numeric_type = T)

    if test
        @test problem_vexity(p) == ConvexVexity()
    end
    handle_problem!(p)
    if test
        @test objective_value(p) === 1.0
        @test evaluate(sum(neg(x))) ≈ 6 atol = atol rtol = rtol
    end
end

@add_problem lp function lp_sumlargest_atom(
    handle_problem!,
    ::Val{test},
    atol,
    rtol,
    ::Type{T},
) where {T,test}
    x = Variable(2)
    p = minimize(
        sumlargest(x, 2) + sumlargest(x, 0),
        x >= [1; 1];
        numeric_type = T,
    )

    if test
        @test problem_vexity(p) == ConvexVexity()
    end
    handle_problem!(p)
    if test
        @test objective_value(p) ≈ 2 atol = atol rtol = rtol
        @test evaluate(sumlargest(x, 2)) ≈ 2 atol = atol rtol = rtol
    end

    x = Variable(4, 4)
    p = minimize(
        sumlargest(x, 3),
        x >= eye(4),
        x[1, 1] >= 1.5,
        x[2, 3] >= 2.1;
        numeric_type = T,
    )

    if test
        @test problem_vexity(p) == ConvexVexity()
    end
    handle_problem!(p)
    if test
        @test objective_value(p) ≈ 4.6 atol = atol rtol = rtol
        @test evaluate(sumlargest(x, 2)) ≈ 3.6 atol = atol rtol = rtol
    end
end

@add_problem lp function lp_sumsmallest_atom(
    handle_problem!,
    ::Val{test},
    atol,
    rtol,
    ::Type{T},
) where {T,test}
    x = Variable(4, 4)
    p = minimize(
        sumlargest(x, 2) + sumsmallest(x, 0),
        sumsmallest(x, 4) >= 1;
        numeric_type = T,
    )

    if test
        @test problem_vexity(p) == ConvexVexity()
    end
    handle_problem!(p)
    if test
        @test objective_value(p) ≈ 0.5 atol = atol rtol = rtol
        @test evaluate(sumsmallest(x, 4)) ≈ 1 atol = atol rtol = rtol
    end

    x = Variable(3, 2)
    p = maximize(
        sumsmallest(x, 3),
        x >= 2,
        x <= 5,
        sumlargest(x, 3) <= 12;
        numeric_type = T,
    )

    if test
        @test problem_vexity(p) == ConvexVexity()
    end
    handle_problem!(p)
    if test
        @test objective_value(p) ≈ 12 atol = atol rtol = rtol
        @test evaluate(sumsmallest(x, 3)) ≈ 12 atol = atol rtol = rtol
    end
end

@add_problem lp function lp_dotsort_atom(
    handle_problem!,
    ::Val{test},
    atol,
    rtol,
    ::Type{T},
) where {T,test}
    x = Variable(4, 1)
    p = minimize(
        dotsort(x, [1, 2, 3, 4]),
        sum(x) >= 7,
        x >= 0,
        x <= 2,
        x[4] <= 1;
        numeric_type = T,
    )

    if test
        @test problem_vexity(p) == ConvexVexity()
    end
    handle_problem!(p)
    if test
        @test objective_value(p) ≈ 19 atol = atol rtol = rtol
        @test vec(evaluate(x)) ≈ [2; 2; 2; 1] atol = atol rtol = rtol
        @test evaluate(dotsort(x, [1, 2, 3, 4])) ≈ 19 atol = atol rtol = rtol
    end

    x = Variable(2, 2)
    p = minimize(
        dotsort(x, [1 2; 3 4]),
        sum(x) >= 7,
        x >= 0,
        x <= 2,
        x[2, 2] <= 1;
        numeric_type = T,
    )

    if test
        @test problem_vexity(p) == ConvexVexity()
    end
    handle_problem!(p)
    if test
        @test objective_value(p) ≈ 19 atol = atol rtol = rtol
        @test evaluate(dotsort(x, [1, 2, 3, 4])) ≈ 19 atol = atol rtol = rtol
    end
end

@add_problem lp function lp_hinge_loss_atom(
    handle_problem!,
    ::Val{test},
    atol,
    rtol,
    ::Type{T},
) where {T,test}
    # TODO: @davidlizeng. We should finish this someday.
end

@add_problem lp function lp_dual_norm_inf_atom(
    handle_problem!,
    ::Val{test},
    atol,
    rtol,
    ::Type{T},
) where {T,test}
    x = Variable(3)
    p = minimize(norm(x, Inf), [-2 <= x, x <= 1]; numeric_type = T)

    if test
        @test problem_vexity(p) == ConvexVexity()
    end
    handle_problem!(p)
    if test
        @test objective_value(p) ≈ 0 atol = atol rtol = rtol
        @test evaluate(norm(x, Inf)) ≈ 0 atol = atol rtol = rtol
        @test norm(p.constraints[1].dual) ≈ 0 atol = atol rtol = rtol
        @test norm(p.constraints[2].dual) ≈ 0 atol = atol rtol = rtol
    end
end

@add_problem lp function lp_dual_norm_1_atom(
    handle_problem!,
    ::Val{test},
    atol,
    rtol,
    ::Type{T},
) where {T,test}
    x = Variable(3)
    p = minimize(norm(x, 1), [-2 <= x, x <= 1]; numeric_type = T)

    if test
        @test problem_vexity(p) == ConvexVexity()
    end
    handle_problem!(p)
    if test
        @test objective_value(p) ≈ 0 atol = atol rtol = rtol
        @test evaluate(norm(x, 1)) ≈ 0 atol = atol rtol = rtol
        @test norm(p.constraints[1].dual) ≈ 0 atol = atol rtol = rtol
        @test norm(p.constraints[2].dual) ≈ 0 atol = atol rtol = rtol
    end
end
