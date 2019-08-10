@testset "Constant variables: $solver" for solver in solvers

    @testset "Issue #166" begin
        # Issue #166
        α = Variable(5)
        fix!(α, ones(5,1))

        # has const vexity, but not at the head
        c = (rand(5,5) * α) * ones(1,5) 

        β = Variable(5)
        β.value = ones(5)

        problem = minimize(sum(c * β), [β >= 0])
        solve!(problem, solver)
        @test problem.optval ≈ evaluate(sum(c * β)) atol=TOL
        @test problem.optval ≈ 0.0 atol=TOL
        @test β.value ≈ zeros(5) atol=TOL
    end

    @testset "Issue #228" begin
        x = Variable(2)
        y = Variable(2)
        fix!(x, [1 1]')
        prob = minimize(y'*(x+[2 2]'), [y>=0])
        solve!(prob, solver)
        @test prob.optval ≈ 0.0 atol = TOL

        prob = minimize(x'*y, [y>=0])
        solve!(prob, solver)
        @test prob.optval ≈ 0.0 atol = TOL
    end

    @testset "Test double `fix!`" begin
        x = Variable()
        y = Variable()
        fix!(x, 1.0)
        prob = minimize(y*x, [y >= x, x >= 0.5])
        solve!(prob, solver)
        @test prob.optval ≈ 1.0 atol = TOL

        fix!(x, 2.0)
        solve!(prob, solver)
        @test prob.optval ≈ 4.0 atol = TOL

        free!(x)
        fix!(y, 1.0)
        solve!(prob, solver)
        @test prob.optval ≈ 0.5 atol = TOL
    end

    @testset "fix! and multiply" begin
        p = Variable()
        fix!(p, 1.0)
        x = Variable(2,2)
        prob = minimize( tr(p*x), [ x >= 1 ])
        solve!(prob, solver)
        @test prob.optval ≈ 2.0 atol = TOL
        @test evaluate( tr(p*x) ) ≈ 2.0 atol = TOL
    end

    @testset "fix! with complex numbers" begin
        x = ComplexVariable()
        fix!(x, 1.0 + im*1.0)
        y = Variable()
        prob = minimize( real(x*y), [ y >= .5, real(x) >= .5, imag(x) >= 0])
        solve!(prob, solver)
        @test prob.optval ≈ .5 atol=TOL
        @test evaluate(real(x*y)) ≈ .5 atol=TOL
        @test evaluate(y) ≈ 0.5 atol=TOL

        free!(x)
        fix!(y)
        solve!(prob, solver)
        @test prob.optval ≈ 0.25 atol=TOL
        @test evaluate(real(x*y)) ≈ 0.25 atol=TOL
        @test real(evaluate(x)) ≈ 0.5 atol=TOL
        @test evaluate(y) ≈ 0.5 atol=TOL
    end

end