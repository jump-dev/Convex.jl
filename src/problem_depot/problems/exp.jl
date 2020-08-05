@add_problem exp function exp_exp_atom(handle_problem!, ::Val{test}, atol, rtol, ::Type{T}) where {T, test}
    y = Variable()
    p = minimize(exp(y), y>=0; numeric_type = T)

    if test
        @test vexity(p) == ConvexVexity()
    end
    handle_problem!(p)
    if test
        @test p.optval ≈ 1 atol=atol rtol=rtol
        @test evaluate(exp(y)) ≈ 1 atol=atol rtol=rtol
    end

    y = Variable()
    p = minimize(exp(y), y>=1; numeric_type = T)

    if test
        @test vexity(p) == ConvexVexity()
    end
    handle_problem!(p)
    if test
        @test p.optval ≈ exp(1) atol=atol rtol=rtol
        @test evaluate(exp(y)) ≈ exp(1) atol=atol rtol=rtol
    end

    y = Variable(5)
    p = minimize(sum(exp(y)), y>=0; numeric_type = T)

    if test
        @test vexity(p) == ConvexVexity()
    end
    handle_problem!(p)
    if test
        @test p.optval ≈ 5 atol=atol rtol=rtol
        @test evaluate(sum(exp(y))) ≈ 5 atol=atol rtol=rtol
    end

    y = Variable(5)
    p = minimize(sum(exp(y)), y>=0; numeric_type = T)

    if test
        @test vexity(p) == ConvexVexity()
    end
    handle_problem!(p)
    if test
        @test p.optval ≈ 5 atol=atol rtol=rtol
    end
end

@add_problem exp function exp_log_atom(handle_problem!, ::Val{test}, atol, rtol, ::Type{T}) where {T, test}
    y = Variable()
    p = maximize(log(y), y<=1; numeric_type = T)

    if test
        @test vexity(p) == ConvexVexity()
    end
    handle_problem!(p)
    if test
        @test p.optval ≈ 0 atol=atol rtol=rtol
    end

    y = Variable()
    p = maximize(log(y), y<=2; numeric_type = T)

    if test
        @test vexity(p) == ConvexVexity()
    end
    handle_problem!(p)
    if test
        @test p.optval ≈ log(2) atol=atol rtol=rtol
    end

    y = Variable()
    p = maximize(log(y), [y<=2, exp(y)<=10]; numeric_type = T)

    if test
        @test vexity(p) == ConvexVexity()
    end
    handle_problem!(p)
    if test
        @test p.optval ≈ log(2) atol=atol rtol=rtol
    end
end

@add_problem exp function exp_log_sum_exp_atom(handle_problem!, ::Val{test}, atol, rtol, ::Type{T}) where {T, test}
    y = Variable(5)
    p = minimize(logsumexp(y), y >= 1; numeric_type = T)

    if test
        @test vexity(p) == ConvexVexity()
    end
    handle_problem!(p)
    if test
        @test p.optval ≈ log(exp(1) * 5) atol=atol rtol=rtol
    end
end

@add_problem exp function exp_logistic_loss_atom(handle_problem!, ::Val{test}, atol, rtol, ::Type{T}) where {T, test}
    y = Variable(5)
    p = minimize(logisticloss(y), y>=1; numeric_type = T)

    if test
        @test vexity(p) == ConvexVexity()
    end
    handle_problem!(p)
    if test
        @test p.optval ≈ log(exp(1) + 1) * 5 atol=atol rtol=rtol
    end
end

@add_problem exp function exp_entropy_atom(handle_problem!, ::Val{test}, atol, rtol, ::Type{T}) where {T, test}
    y = Variable(5, Positive())
    p = maximize(entropy(y), sum(y)<=1; numeric_type = T)

    if test
        @test vexity(p) == ConvexVexity()
    end
    handle_problem!(p)
    if test
        @test p.optval ≈ -(log(1 / 5)) atol=atol rtol=rtol
    end
end

@add_problem exp function exp_relative_entropy_atom(handle_problem!, ::Val{test}, atol, rtol, ::Type{T}) where {T, test}
    x = Variable(1)
    y = Variable(1)
    # x log (x/y)
    p = minimize(relative_entropy(x,y), y==1, x >= 2; numeric_type = T)

    if test
        @test vexity(p) == ConvexVexity()
    end
    handle_problem!(p)
    if test
        @test p.optval ≈ 2 * log(2) atol=atol rtol=rtol
    end
end

@add_problem exp function exp_log_perspective_atom(handle_problem!, ::Val{test}, atol, rtol, ::Type{T}) where {T, test}
    x = Variable(1)
    y = Variable(1)
    # y log (x/y)
    p = maximize(log_perspective(x,y), y==5, x <= 10; numeric_type = T)

    if test
        @test vexity(p) == ConvexVexity()
    end
    handle_problem!(p)
    if test
        @test p.optval ≈ 5 * log(2) atol=atol rtol=rtol
    end
end
