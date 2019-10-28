@add_problem mip function mip_lp_fallback_interface(handle_problem!, ::Val{test}, atol, rtol, ::Type{T}) where {T, test}
    x = Variable()
    p = minimize(x, x>=4.3)
    if test
        @test vexity(p) == AffineVexity()
    end
    handle_problem!(p)
    if test
        @test p.optval ≈ 4.3 atol=atol rtol=rtol
    end

    x = Variable(2)
    p = minimize(norm(x,1), x[1]>=4.3)
    if test
        @test vexity(p) == ConvexVexity()
    end
    handle_problem!(p)
    if test
        @test p.optval ≈ 4.3 atol=atol rtol=rtol
    end
end

@add_problem mip function mip_integer_variables(handle_problem!, ::Val{test}, atol, rtol, ::Type{T}) where {T, test}
    x = Variable(:Int)
    p = minimize(x, x>=4.3)
    if test
        @test vexity(p) == AffineVexity()
    end
    handle_problem!(p)
    if test
        @test p.optval ≈ 5 atol=atol rtol=rtol
    end

    x = Variable(2, :Int)
    p = minimize(sum(x), x>=4.3)
    if test
        @test vexity(p) == AffineVexity()
    end
    handle_problem!(p)
    if test
        @test p.optval ≈ 10 atol=atol rtol=rtol
    end

    x = Variable(:Int)
    y = Variable()
    p = minimize(sum(x + y), x>=4.3, y>=7)
    if test
        @test vexity(p) == AffineVexity()
    end
    handle_problem!(p)
    if test
        @test p.optval ≈ 12 atol=atol rtol=rtol
    end

    x = Variable(2, :Int)
    p = minimize(norm(x, 1), x[1]>=4.3)
    if test
        @test vexity(p) == ConvexVexity()
    end
    handle_problem!(p)
    if test
        @test p.optval ≈ 5 atol=atol rtol=rtol
    end

    x = Variable(2, :Int)
    p = minimize(sum(x), x[1]>=4.3, x>=0)
    if test
        @test vexity(p) == AffineVexity()
    end
    handle_problem!(p)
    if test
        @test p.optval ≈ 5 atol=atol rtol=rtol
    end

    x = Variable(2, :Int)
    p = minimize(sum(x), x>=.5)
    if test
        @test vexity(p) == AffineVexity()
    end
    handle_problem!(p)
    if test
        @test p.optval ≈ 2 atol=atol rtol=rtol
    end
end

@add_problem mip function mip_binary_variables(handle_problem!, ::Val{test}, atol, rtol, ::Type{T}) where {T, test}
    x = Variable(2, :Bin)
    p = minimize(sum(x), x>=.5)
    if test
        @test vexity(p) == AffineVexity()
    end
    handle_problem!(p)
    if test
        @test p.optval ≈ 2 atol=atol rtol=rtol
    end

    x = Variable(2, :Bin)
    p = minimize(sum(x), x[1]>=.5, x>=0)
    if test
        @test vexity(p) == AffineVexity()
    end
    handle_problem!(p)
    if test
        @test p.optval ≈ 1 atol=atol rtol=rtol
    end
end
