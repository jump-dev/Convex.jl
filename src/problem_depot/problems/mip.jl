@add_problem mip function mip_lp_fallback_interface(
    handle_problem!,
    ::Val{test},
    atol,
    rtol,
    ::Type{T},
) where {T,test}
    x = Variable()
    p = minimize(x, x >= 4.3; numeric_type = T)

    if test
        @test problem_vexity(p) == AffineVexity()
    end
    handle_problem!(p)
    if test
        @test objective_value(p) ≈ 4.3 atol = atol rtol = rtol
    end

    x = Variable(2)
    p = minimize(norm(x, 1), x[1] >= 4.3; numeric_type = T)

    if test
        @test problem_vexity(p) == ConvexVexity()
    end
    handle_problem!(p)
    if test
        @test objective_value(p) ≈ 4.3 atol = atol rtol = rtol
    end
end

@add_problem mip function mip_integer_variables(
    handle_problem!,
    ::Val{test},
    atol,
    rtol,
    ::Type{T},
) where {T,test}
    x = Variable(:Int)
    p = minimize(x, x >= 4.3; numeric_type = T)

    if test
        @test problem_vexity(p) == AffineVexity()
    end
    handle_problem!(p)
    if test
        @test objective_value(p) ≈ 5 atol = atol rtol = rtol
    end

    x = Variable(2, IntVar)
    p = minimize(sum(x), x >= 4.3; numeric_type = T)

    if test
        @test problem_vexity(p) == AffineVexity()
    end
    handle_problem!(p)
    if test
        @test objective_value(p) ≈ 10 atol = atol rtol = rtol
    end

    x = Variable(:Int)
    y = Variable()
    p = minimize(sum(x + y), x >= 4.3, y >= 7; numeric_type = T)

    if test
        @test problem_vexity(p) == AffineVexity()
    end
    handle_problem!(p)
    if test
        @test objective_value(p) ≈ 12 atol = atol rtol = rtol
    end

    x = Variable(2, IntVar)
    p = minimize(norm(x, 1), x[1] >= 4.3; numeric_type = T)

    if test
        @test problem_vexity(p) == ConvexVexity()
    end
    handle_problem!(p)
    if test
        @test objective_value(p) ≈ 5 atol = atol rtol = rtol
    end

    x = Variable(2, IntVar)
    p = minimize(sum(x), x[1] >= 4.3, x >= 0; numeric_type = T)

    if test
        @test problem_vexity(p) == AffineVexity()
    end
    handle_problem!(p)
    if test
        @test objective_value(p) ≈ 5 atol = atol rtol = rtol
    end

    x = Variable(2, IntVar)
    p = minimize(sum(x), x >= 0.5; numeric_type = T)

    if test
        @test problem_vexity(p) == AffineVexity()
    end
    handle_problem!(p)
    if test
        @test objective_value(p) ≈ 2 atol = atol rtol = rtol
    end
end

@add_problem mip function mip_binary_variables(
    handle_problem!,
    ::Val{test},
    atol,
    rtol,
    ::Type{T},
) where {T,test}
    x = Variable(2, BinVar)
    p = minimize(sum(x), x >= 0.5; numeric_type = T)

    if test
        @test problem_vexity(p) == AffineVexity()
    end
    handle_problem!(p)
    if test
        @test objective_value(p) ≈ 2 atol = atol rtol = rtol
    end

    x = Variable(2, BinVar)
    p = minimize(sum(x), x[1] >= 0.5, x >= 0; numeric_type = T)

    if test
        @test problem_vexity(p) == AffineVexity()
    end
    handle_problem!(p)
    if test
        @test objective_value(p) ≈ 1 atol = atol rtol = rtol
    end
end
