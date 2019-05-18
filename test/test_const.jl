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


end