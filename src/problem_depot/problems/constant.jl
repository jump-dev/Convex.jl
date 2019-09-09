

@add_problem constant function constant_Issue_166(handle_problem!, ::Val{test}, atol, rtol, ::Type{T}) where {T, test}
    # Issue #166
    α = Variable(5)
    fix!(α, ones(5,1))

    # has const vexity, but not at the head
    c = (rand(5,5) * α) * ones(1,5) 

    β = Variable(5)
    β.value = ones(5)

    problem = minimize(sum(c * β), [β >= 0])
    handle_problem!(problem)
    if test
        @test problem.optval ≈ evaluate(sum(c * β)) atol=atol rtol=rtol
        @test problem.optval ≈ 0.0 atol=atol rtol=rtol
        @test β.value ≈ zeros(5) atol=atol rtol=rtol
    end
end

@add_problem constant function constant_Issue_228(handle_problem!, ::Val{test}, atol, rtol, ::Type{T}) where {T, test}
    x = Variable(2)
    y = Variable(2)
    fix!(x, [1 1]')
    prob = minimize(y'*(x+[2 2]'), [y>=0])
    handle_problem!(prob)
    if test
        @test prob.optval ≈ 0.0 atol=atol rtol=rtol
    end

    prob = minimize(x'*y, [y>=0])
    handle_problem!(prob)
    if test
        @test prob.optval ≈ 0.0 atol=atol rtol=rtol
    end
end

@add_problem constant function constant_Test_double_fix!(handle_problem!, ::Val{test}, atol, rtol, ::Type{T}) where {T, test}
    x = Variable()
    y = Variable()
    fix!(x, 1.0)
    prob = minimize(y*x, [y >= x, x >= 0.5])
    handle_problem!(prob)
    if test
        @test prob.optval ≈ 1.0 atol=atol rtol=rtol
    end

    fix!(x, 2.0)
    handle_problem!(prob)
    if test
        @test prob.optval ≈ 4.0 atol=atol rtol=rtol
    end

    free!(x)
    fix!(y, 1.0)
    handle_problem!(prob)
    if test
        @test prob.optval ≈ 0.5 atol=atol rtol=rtol
    end
end

@add_problem constant function constant_fix!_and_multiply(handle_problem!, ::Val{test}, atol, rtol, ::Type{T}) where {T, test}
    p = Variable()
    fix!(p, 1.0)
    x = Variable(2,2)
    prob = minimize( tr(p*x), [ x >= 1 ])
    handle_problem!(prob)
    if test
        @test prob.optval ≈ 2.0 atol=atol rtol=rtol
        @test evaluate( tr(p*x) ) ≈ 2.0 atol=atol rtol=rtol
    end
end

@add_problem constant function constant_fix!_with_complex_numbers(handle_problem!, ::Val{test}, atol, rtol, ::Type{T}) where {T, test}
    x = ComplexVariable()
    fix!(x, 1.0 + im*1.0)
    y = Variable()
    prob = minimize( real(x*y), [ y >= .5, real(x) >= .5, imag(x) >= 0])
    handle_problem!(prob)
    if test
        @test prob.optval ≈ .5 atol=atol rtol=rtol
        @test evaluate(real(x*y)) ≈ .5 atol=atol rtol=rtol
        @test evaluate(y) ≈ 0.5 atol=atol rtol=rtol
    end

    free!(x)
    fix!(y)
    handle_problem!(prob)
    if test
        @test prob.optval ≈ 0.25 atol=atol rtol=rtol
        @test evaluate(real(x*y)) ≈ 0.25 atol=atol rtol=rtol
        @test real(evaluate(x)) ≈ 0.5 atol=atol rtol=rtol
        @test evaluate(y) ≈ 0.5 atol=atol rtol=rtol
    end

    if test
        @test_throws DimensionMismatch fix!(x, rand(2))
        @test_throws DimensionMismatch fix!(x, rand(2,2))
    end
end

@add_problem constant function constant_fix!_with_vectors(handle_problem!, ::Val{test}, atol, rtol, ::Type{T}) where {T, test}
    x = ComplexVariable(5)
    fix!(x, ones(5) + im*ones(5))
    y = Variable()
    prob = minimize( real(y*sum(x)), [ y >= .5, real(x) >= .5, imag(x) >= 0])
    handle_problem!(prob)
    if test
        @test prob.optval ≈ 2.5 atol=atol rtol=rtol
        @test evaluate(real(y*sum(x))) ≈ 2.5 atol=atol rtol=rtol
        @test evaluate(y) ≈ 0.5 atol=atol rtol=rtol
    end

    free!(x)
    fix!(y)
    handle_problem!(prob)
    if test
        @test prob.optval ≈ 1.25 atol=atol rtol=rtol
        @test evaluate(real(y*sum(x))) ≈ 1.25 atol=atol rtol=rtol
        @test real(evaluate(x)) ≈ 0.5*ones(5) atol=atol rtol=rtol
        @test evaluate(y) ≈ 0.5 atol=atol rtol=rtol
    end

    if test
        @test_throws DimensionMismatch fix!(x, rand(5,5))
        @test_throws DimensionMismatch fix!(x, rand(4))
    end

end
