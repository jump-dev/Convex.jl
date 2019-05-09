@testset "Affine Atoms: $solver" for solver in solvers
    @testset "negate atom" begin
        x = Variable()
        p = minimize(-x, [x <= 0])
        @test vexity(p) == AffineVexity()
        solve!(p, solver)
        @test p.optval ≈ 0 atol=TOL
        @test evaluate(-x) ≈ 0 atol=TOL
    end

    @testset "kron atom" begin
        x = ComplexVariable(3, 3)
        y = [1.0 2.0; 3.0 4.0]
        @test size(kron(x, y)) == (6, 6)
        @test size(kron(y, x)) == (6, 6)
    end

    @testset "multiply atom" begin
        x = Variable(1)
        p = minimize(2.0 * x, [x >= 2, x <= 4])
        @test vexity(p) == AffineVexity()
        solve!(p, solver)
        @test p.optval ≈ 4 atol=TOL
        @test (evaluate(2.0x))[1] ≈ 4 atol=TOL

        x = Variable(2)
        A = 1.5 * eye(2)
        p = minimize([2 2] * x, [A * x >= [1.1; 1.1]])
        @test vexity(p) == AffineVexity()
        solve!(p, solver)
        @test p.optval ≈ 2.93333 atol=TOL
        @test (evaluate([2 2] * x))[1] ≈ 2.93333 atol=TOL
        @test vec(evaluate(A * x)) ≈ [1.1; 1.1] atol=TOL

        y = Variable(1)
        x = Variable(3)
        z = [1.0, 2.0, 3.0] * y
        k = -y * [1.0, 2.0, 3.0]
        c = [y <= 3.0, y >= 0.0, x >= ones(3), k <= x, x <= z]
        o = 3 * y
        p = Problem(:minimize, o, c)
        @test vexity(p) == AffineVexity()
        solve!(p, solver)
        @test p.optval ≈ 3 atol=TOL

        p = Problem(:minimize, o, c...)
        @test vexity(p) == AffineVexity()
        solve!(p, solver)
        @test p.optval ≈ 3 atol=TOL

        # Check #274
        x = ComplexVariable(2,2)
        p = minimize( real( [1.0im, 0.0]' * x * [1.0im, 0.0] ), [ x == [1.0 0.0; 0.0 1.0] ])
        solve!(p, solver)
        @test p.optval ≈ 1.0
    end

    @testset "dot atom" begin
        x = Variable(2)
        p = minimize(dot([2.0; 2.0], x), x >= [1.1; 1.1])
        @test vexity(p) == AffineVexity()
        solve!(p, solver)
        @test p.optval ≈ 4.4 atol=TOL
        @test (evaluate(dot([2.0; 2.0], x)))[1] ≈ 4.4 atol=TOL
    end

    @testset "dot atom for matrix variables" begin
        x = Variable(2,2)
        p = minimize(dot(fill(2.0, (2,2)), x), x >= 1.1)
        @test vexity(p) == AffineVexity()
        solve!(p, solver)
        @test p.optval ≈ 8.8 atol=TOL
        @test (evaluate(dot(fill(2.0, (2, 2)), x)))[1] ≈ 8.8 atol=TOL
    end

    @testset "add atom" begin
        x = Variable(1)
        y = Variable(1)
        p = minimize(x + y, [x >= 3, y >= 2])
        @test vexity(p) == AffineVexity()
        solve!(p, solver)
        @test p.optval ≈ 5 atol=TOL
        @test evaluate(x + y) ≈ 5 atol=TOL

        x = Variable(1)
        p = minimize(x, [eye(2) + x >= eye(2)])
        @test vexity(p) == AffineVexity()
        solve!(p, solver)
        @test p.optval ≈ 0 atol=TOL
        @test evaluate(eye(2) + x) ≈ eye(2) atol=TOL

        y = Variable()
        p = minimize(y - 5, y >= -1)
        @test vexity(p) == AffineVexity()
        solve!(p, solver)
        @test p.optval ≈ -6 atol=TOL
        @test evaluate(y - 5) ≈ -6 atol=TOL
    end

    @testset "transpose atom" begin
        x = Variable(2)
        c = ones(2, 1)
        p = minimize(x' * c, x >= 1)
        @test vexity(p) == AffineVexity()
        solve!(p, solver)
        @test p.optval ≈ 2 atol=TOL
        @test (evaluate(x' * c))[1] ≈ 2 atol=TOL

        X = Variable(2, 2)
        c = ones(2, 1)
        p = minimize(c' * X' * c, [X >= ones(2, 2)])
        @test vexity(p) == AffineVexity()
        solve!(p, solver)
        @test p.optval ≈ 4 atol=TOL
        @test (evaluate(c' * X' * c))[1] ≈ 4 atol=TOL

        rows = 2
        cols = 3
        r = rand(rows, cols)
        r_2 = rand(cols, rows)
        x = Variable(rows, cols)
        c = ones(1, cols)
        d = ones(rows, 1)
        p = minimize(c * x' * d + d' * x * c' + (c * x''''' * d)',
                    [x' >= r_2, x >= r, x''' >= r_2, x'' >= r])
        @test vexity(p) == AffineVexity()
        solve!(p, solver)
        s = sum(max.(r, r_2')) * 3
        @test p.optval ≈ s atol=TOL
        @test (evaluate(c * x' * d + d' * x * c' + (c * ((((x')')')')' * d)'))[1] ≈ s atol=TOL
    end

    @testset "index atom" begin
        x = Variable(2)
        p = minimize(x[1] + x[2], [x >= 1])
        @test vexity(p) == AffineVexity()
        solve!(p, solver)
        @test p.optval ≈ 2 atol=TOL
        @test (evaluate(x[1] + x[2]))[1] ≈ 2 atol=TOL

        x = Variable(3)
        I = [true true false]
        p = minimize(sum(x[I]), [x >= 1])
        @test vexity(p) == AffineVexity()
        solve!(p, solver)
        @test p.optval ≈ 2 atol=TOL
        @test (evaluate(sum(x[I])))[1] ≈ 2 atol=TOL

        rows = 6
        cols = 8
        n = 2
        X = Variable(rows, cols)
        A = randn(rows, cols)
        c = rand(1, n)
        p = minimize(c * X[1:n, 5:5+n-1]' * c', X >= A)
        @test vexity(p) == AffineVexity()
        solve!(p, solver)
        s = c * A[1:n, 5:5+n-1]' * c'
        @test p.optval ≈ s[1] atol=TOL
        @test evaluate(c * (X[1:n, 5:(5 + n) - 1])' * c') ≈ s atol=TOL
    end

    @testset "sum atom" begin
        x = Variable(2,2)
        p = minimize(sum(x), x>=1)
        @test vexity(p) == AffineVexity()
        solve!(p, solver)
        @test p.optval ≈ 4 atol=TOL
        @test evaluate(sum(x)) ≈ 4 atol=TOL

        x = Variable(2,2)
        p = minimize(sum(x) - 2*x[1,1], x>=1, x[1,1]<=2)
        @test vexity(p) == AffineVexity()
        solve!(p, solver)
        @test p.optval ≈ 1 atol=TOL
        @test (evaluate(sum(x) - 2 * x[1, 1]))[1] ≈ 1 atol=TOL

        x = Variable(10)
        a = rand(10, 1)
        p = maximize(sum(x[2:6]), x <= a)
        @test vexity(p) == AffineVexity()
        solve!(p, solver)
        @test p.optval ≈ sum(a[2:6]) atol=TOL
        @test evaluate(sum(x[2:6])) ≈ sum(a[2:6]) atol=TOL
    end

    @testset "diag atom" begin
        x = Variable(2,2)
        p = minimize(sum(diag(x,1)), x >= 1)
        @test vexity(p) == AffineVexity()
        solve!(p, solver)
        @test p.optval ≈ 1 atol=TOL
        @test evaluate(sum(diag(x, 1))) ≈ 1 atol=TOL

        x = Variable(4, 4)
        p = minimize(sum(diag(x)), x >= 2)
        @test vexity(p) == AffineVexity()
        solve!(p, solver)
        @test p.optval ≈ 8 atol=TOL
        @test evaluate(sum(diag(x))) ≈ 8 atol=TOL
    end

    @testset "trace atom" begin
        x = Variable(2,2)
        p = minimize(tr(x), x >= 1)
        @test vexity(p) == AffineVexity()
        solve!(p, solver)
        @test p.optval ≈ 2 atol=TOL
        @test evaluate(tr(x)) ≈ 2 atol=TOL
    end

    @testset "dot multiply atom" begin
        x = Variable(3)
        p = maximize(sum(dot(*)(x,[1,2,3])), x<=1)
        @test vexity(p) == AffineVexity()
        solve!(p, solver)
        @test p.optval ≈ 6 atol=TOL
        @test evaluate(sum((dot(*))(x, [1, 2, 3]))) ≈ 6 atol=TOL

        x = Variable(3, 3)
        p = maximize(sum(dot(*)(x,eye(3))), x<=1)
        @test vexity(p) == AffineVexity()
        solve!(p, solver)
        @test p.optval ≈ 3 atol=TOL
        @test evaluate(sum((dot(*))(x, eye(3)))) ≈ 3 atol=TOL

        x = Variable(5, 5)
        p = minimize(x[1, 1], dot(*)(3,x) >= 3)
        @test vexity(p) == AffineVexity()
        solve!(p, solver)
        @test p.optval ≈ 1 atol=TOL
        @test (evaluate(x[1, 1]))[1] ≈ 1 atol=TOL

        x = Variable(3,1)
        p = minimize(sum(dot(*)(ones(3,3), x)), x>=1)
        @test vexity(p) == AffineVexity()
        solve!(p, solver)
        @test p.optval ≈ 9 atol=TOL
        @test (evaluate(x[1, 1]))[1] ≈ 1 atol=TOL

        x = Variable(1,3)
        p = minimize(sum(dot(*)(ones(3,3), x)), x>=1)
        @test vexity(p) == AffineVexity()
        solve!(p, solver)
        @test p.optval ≈ 9 atol=TOL
        @test (evaluate(x[1, 1]))[1] ≈ 1 atol=TOL

        x = Variable(1, 3, Positive())
        p = maximize(sum(dot(/)(x,[1 2 3])), x<=1)
        @test vexity(p) == AffineVexity()
        solve!(p, solver)
        @test p.optval ≈ 11 / 6 atol=TOL
        @test evaluate(sum((dot(/))(x, [1 2 3]))) ≈ 11 / 6 atol=TOL

        # Broadcast fusion works
        x = Variable(5, 5)
        a = 2.0 .* x .* ones(Int, 5)
        @test a isa DotMultiplyAtom
    end

    @testset "reshape atom" begin
        A = rand(2, 3)
        X = Variable(3, 2)
        c = rand()
        p = minimize(sum(reshape(X, 2, 3) + A), X >= c)
        @test vexity(p) == AffineVexity()
        solve!(p, solver)
        @test p.optval ≈ sum(A .+ c) atol=TOL
        @test evaluate(sum(reshape(X, 2, 3) + A)) ≈ sum(A .+ c) atol=TOL

        b = rand(6)
        p = minimize(sum(vec(X) + b), X >= c)
        @test vexity(p) == AffineVexity()
        solve!(p, solver)
        @test p.optval ≈ sum(b .+ c) atol=TOL
        @test evaluate(sum(vec(X) + b)) ≈ sum(b .+ c) atol=TOL

        x = Variable(4, 4)
        c = ones(16, 1)
        reshaped = reshape(x, 16, 1)
        a = collect(1:16)
        p = minimize(c' * reshaped, reshaped >= a)
        @test vexity(p) == AffineVexity()
        solve!(p, solver)
        # TODO: why is accuracy lower here?
        @test p.optval ≈ 136 atol=10TOL
        @test (evaluate(c' * reshaped))[1] ≈ 136 atol=10TOL
    end

    @testset "hcat atom" begin
        x = Variable(4, 4)
        y = Variable(4, 6)
        p = maximize(sum(x) + sum([y fill(4.0, 4)]), [x y fill(2.0, (4, 2))] <= 2)
        @test vexity(p) == AffineVexity()
        solve!(p, solver)
        @test p.optval ≈ 96 atol=TOL
        @test evaluate(sum(x) + sum([y fill(4.0, 4)])) ≈ 96 atol=TOL
        @test evaluate([x y fill(2.0, (4, 2))]) ≈ fill(2.0, (4, 12)) atol=TOL
    end

    @testset "vcat atom" begin
        x = Variable(4, 4)
        y = Variable(4, 6)

        # TODO: fix dimension mismatch [y 4*eye(4); x -ones(4, 6)]
        p = maximize(sum(x) + sum([y 4*eye(4); x -ones(4, 6)]), [x;y'] <= 2)
        @test vexity(p) == AffineVexity()
        solve!(p, solver)
        # TODO: why is accuracy lower here?
        @test p.optval ≈ 104 atol=10TOL
        @test evaluate(sum(x) + sum([y 4 * eye(4); x -(ones(4, 6))])) ≈ 104 atol=10TOL
        @test evaluate([x; y']) ≈ 2 * ones(10, 4) atol=TOL
    end

    @testset "Diagonal atom" begin
        x = Variable(2, 2)
        @test_throws ArgumentError Diagonal(x)

        x = Variable(4)
        p = minimize(sum(Diagonal(x)), x == [1, 2, 3, 4])
        @test vexity(p) == AffineVexity()
        solve!(p, solver)
        @test p.optval ≈ 10 atol=TOL
        @test all(abs.(evaluate(Diagonal(x)) - Diagonal([1, 2, 3, 4])) .<= TOL)

        x = Variable(3)
        c = [1; 2; 3]
        p = minimize(c' * Diagonal(x) * c, x >= 1, sum(x) == 10)
        @test vexity(p) == AffineVexity()
        solve!(p, solver)
        @test p.optval ≈ 21 atol=TOL

        x = Variable(3)
        p = minimize(sum(x), x >= 1, Diagonal(x)[1, 2] == 1)
        @test solve!(p, solver) === nothing
        @test p.status != :Optimal
    end

    @testset "conv atom" begin
        x = Variable(3)
        h = [1, -1]
        p = minimize(sum(conv(h, x)) + sum(x), x >= 1, x <= 2)
        @test vexity(p) == AffineVexity()
        solve!(p, solver)
        @test p.optval ≈ 3 atol=TOL
        @test evaluate(sum(conv(h, x))) ≈ 0 atol=TOL

        x = Variable(3)
        h = [1, -1]
        p = minimize(sum(conv(x, h)) + sum(x), x >= 1, x <= 2)
        @test vexity(p) == AffineVexity()
        solve!(p, solver)
        @test p.optval ≈ 3 atol=TOL
        @test evaluate(sum(conv(h, x))) ≈ 0 atol=TOL


        # test #288
        x = rand(5) + im*rand(5);
        y = Constant(rand(5));
        @test conv(x,y) isa AbstractExpr

        x =  Variable(3)
        y = im * x
        h = [1, -1]
        p = minimize(imag(sum(conv(y, h)) + sum(y)), x >= 1, x <= 2)
        @test vexity(p) == AffineVexity()
        solve!(p, solver)
        @test p.optval ≈ 3 atol=TOL
        @test evaluate(sum(conv(h, x))) ≈ 0 atol=TOL

    end

    @testset "satisfy problems" begin
        x = Variable()
        p = satisfy(x >= 0)
        add_constraints!(p, x >= 1)
        add_constraints!(p, [x >= -1, x <= 4])
        solve!(p, solver)
        @test p.status == :Optimal

        p = satisfy([x >= 0, x >= 1, x <= 2])
        solve!(p, solver)
        @test p.status == :Optimal

        p = maximize(1, [x >= 1, x <= 2])
        solve!(p, solver)
        @test p.status == :Optimal

        constr = x >= 0
        constr += x >= 1
        constr += x <= 10
        constr2 = x >= 0
        constr2 += [x >= 2, x <= 3] + constr
        p = satisfy(constr)
        solve!(p, solver)
        @test p.status == :Optimal
    end

    @testset "dual" begin
        x = Variable()
        p = minimize(x, x >= 0)
        solve!(p, solver)
        if p.solution.has_dual
            println("Solution object has dual value, checking for dual correctness.")
            @test p.constraints[1].dual ≈ 1 atol=TOL
        end

        x = Variable()
        p = maximize(x, x <= 0)
        solve!(p, solver)
        if p.solution.has_dual
            println("Solution object has dual value, checking for dual correctness.")
            @test p.constraints[1].dual ≈ 1 atol=TOL
        end

        x = Variable()
        p = minimize(x, x >= 0, x == 2)
        solve!(p, solver)
        if p.solution.has_dual
            println("Solution object has dual value, checking for dual correctness.")
            @test p.constraints[1].dual ≈ 0 atol=TOL
            @test abs.(p.constraints[2].dual) ≈ 1 atol=TOL
        end

        x = Variable(2)
        A = 1.5 * eye(2)
        p = minimize(dot([2.0; 2.0], x), [A * x >= [1.1; 1.1]])
        solve!(p, solver)
        if p.solution.has_dual
            println("Solution object has dual value, checking for dual correctness.")
            dual = [4/3; 4/3]
            @test all(abs.(p.constraints[1].dual - dual) .<= TOL)
        end
    end
end
