using Convex
using Convex: DotMultiplyAtom
using Test
using ECOS
using SCS
using GLPKMathProgInterface
using Random

import LinearAlgebra.eigen
import LinearAlgebra.I
import LinearAlgebra.opnorm
import Random.shuffle
import Statistics.mean

TOL = 1e-3
eye(n) = Matrix(1.0I, n, n)

# Seed random number stream to improve test reliability
Random.seed!(2)

solvers = Any[]

push!(solvers, ECOSSolver(verbose=0))
push!(solvers, GLPKSolverMIP())
push!(solvers, SCSSolver(verbose=0, eps=1e-6))

if isinstalled("Gurobi")
    using Gurobi
    push!(solvers, GurobiSolver(OutputFlag=0))
end

if isinstalled("Mosek")
    using Mosek
    push!(solvers, MosekSolver(LOG=0))
end

@testset "$solver" for solver in solvers
    @testset "utilities" begin
        x = Variable(2,3)
        @test length(x) == 6
        @test size(x) == (2, 3)
        @test size(x, 1) == 2
        @test size(x, 2) == 3

        x = Variable(3)
        @test length(x) == 3
        @test size(x) == (3, 1)

        x = Variable()
        @test length(x) == 1
        @test size(x) == (1, 1)
    end

    @testset "Affine Atoms" begin

        @testset "negate atom" begin
            x = Variable()
            p = minimize(-x, [x <= 0])
            @test vexity(p) == AffineVexity()
            solve!(p, solver)
            @test p.optval ≈ 0 atol=TOL
            @test evaluate(-x) ≈ 0 atol=TOL
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

    @testset "LP Atoms" begin

        @testset "abs atom" begin
            x = Variable()
            p = minimize(abs(x), x<=-1)
            @test vexity(p) == ConvexVexity()
            solve!(p, solver)
            @test p.optval ≈ 1 atol=TOL
            @test evaluate(abs(x)) ≈ 1 atol=TOL

            x = Variable(2,2)
            p = minimize(sum(abs(x)), x[2,2]>=1, x[1,1]>=1, x>=0)
            @test vexity(p) == ConvexVexity()
            solve!(p, solver)
            @test p.optval ≈ 2 atol=TOL
            @test evaluate(sum(abs(x))) ≈ 2 atol=TOL
        end

        @testset "maximum atom" begin
            x = Variable(10)
            a = shuffle(collect(0.1:0.1:1.0))
            p = minimize(maximum(x), x >= a)
            @test vexity(p) == ConvexVexity()
            solve!(p, solver)
            @test p.optval ≈ maximum(a) atol=TOL
            @test evaluate(maximum(x)) ≈ maximum(a) atol=TOL
        end

        @testset "minimum atom" begin
            x = Variable(1)
            a = reshape(shuffle(collect(0.01:0.01:1.0)), (10, 10))
            p = maximize(minimum(x), x <= a)
            @test vexity(p) == ConvexVexity()
            solve!(p, solver)
            @test p.optval ≈ minimum(a) atol=TOL
            @test evaluate(minimum(x)) ≈ minimum(a) atol=TOL

            x = Variable(4, 4)
            y = Variable(4, 6)
            z = Variable(1)
            c = ones(4, 1)
            d = fill(2.0, (6, 1))
            constraints = [[x y] <= 2, z <= 0, z <= x, 2z >= -1]
            objective = sum(x + z) + minimum(y) + c' * y * d
            p = maximize(objective, constraints)
            @test vexity(p) == ConvexVexity()
            solve!(p, solver)
            @test p.optval ≈ 130 atol=TOL
            @test (evaluate(objective))[1] ≈ 130 atol=TOL
        end

        @testset "max atom" begin
            x = Variable(10, 10)
            y = Variable(10, 10)
            a = reshape(shuffle(collect(0.01:0.01:1.0)), (10, 10))
            b = reshape(shuffle(collect(0.01:0.01:1.0)), (10, 10))
            p = minimize(maximum(max(x, y)), [x >= a, y >= b])
            @test vexity(p) == ConvexVexity()
            solve!(p, solver)
            max_a = maximum(a)
            max_b = maximum(b)
            @test p.optval ≈ max(max_a, max_b) atol=10TOL
            @test evaluate(maximum(max(x, y))) ≈ max(max_a, max_b) atol=10TOL
        end

        @testset "min atom" begin
            x = Variable(10, 10)
            y = Variable(10, 10)
            a = reshape(shuffle(collect(0.01:0.01:1.0)), (10, 10))
            b = reshape(shuffle(collect(0.01:0.01:1.0)), (10, 10))
            p = maximize(minimum(min(x, y)), [x <= a, y <= b])
            @test vexity(p) == ConvexVexity()
            solve!(p, solver)
            min_a = minimum(a)
            min_b = minimum(b)
            @test p.optval ≈ min(min_a, min_b) atol=10TOL
            @test evaluate(minimum(min(x, y))) ≈ min(min_a, min_b) atol=10TOL
        end

        @testset "pos atom" begin
            x = Variable(3)
            a = [-2; 1; 2]
            p = minimize(sum(pos(x)), [x >= a, x <= 2])
            @test vexity(p) == ConvexVexity()
            solve!(p, solver)
            @test p.optval ≈ 3 atol=TOL
            @test evaluate(sum(pos(x))) ≈ 3 atol=TOL
        end

        @testset "neg atom" begin
            x = Variable(3)
            p = minimize(1, [x >= -2, x <= -2, neg(x) >= -3])
            @test vexity(p) == ConvexVexity()
            solve!(p, solver)
            @test p.optval ≈ 1 atol=TOL
            @test evaluate(sum(neg(x))) ≈ -6 atol=TOL
        end

        @testset "sumlargest atom" begin
            x = Variable(2)
            p = minimize(sumlargest(x, 2), x >= [1; 1])
            @test vexity(p) == ConvexVexity()
            solve!(p, solver)
            @test p.optval ≈ 2 atol=TOL
            @test evaluate(sumlargest(x, 2)) ≈ 2 atol=TOL

            x = Variable(4, 4)
            p = minimize(sumlargest(x, 3), x >= eye(4), x[1, 1] >= 1.5, x[2, 3] >= 2.1)
            @test vexity(p) == ConvexVexity()
            solve!(p, solver)
            @test p.optval ≈ 4.6 atol=TOL
            @test evaluate(sumlargest(x, 2)) ≈ 3.6 atol=TOL
        end

        @testset "sumsmallest atom" begin
            x = Variable(4, 4)
            p = minimize(sumlargest(x, 2), sumsmallest(x, 4) >= 1)
            @test vexity(p) == ConvexVexity()
            solve!(p, solver)
            @test p.optval ≈ 0.5 atol=TOL
            @test evaluate(sumsmallest(x, 4)) ≈ 1 atol=TOL

            x = Variable(3, 2)
            p = maximize(sumsmallest(x, 3), x >= 2, x <= 5, sumlargest(x, 3) <= 12)
            @test vexity(p) == ConvexVexity()
            solve!(p, solver)
            @test p.optval ≈ 12 atol=TOL
            @test evaluate(sumsmallest(x, 3)) ≈ 12 atol=TOL
        end

        @testset "dotsort atom" begin
            x = Variable(4, 1)
            p = minimize(dotsort(x, [1, 2, 3, 4]), sum(x) >= 7, x >= 0, x <= 2, x[4] <= 1)
            @test vexity(p) == ConvexVexity()
            solve!(p, solver)
            @test p.optval ≈ 19 atol=TOL
            @test vec(x.value) ≈ [2; 2; 2; 1] atol=TOL
            @test evaluate(dotsort(x, [1, 2, 3, 4])) ≈ 19 atol=TOL

            x = Variable(2, 2)
            p = minimize(dotsort(x, [1 2; 3 4]), sum(x) >= 7, x >= 0, x <= 2, x[2, 2] <= 1)
            @test vexity(p) == ConvexVexity()
            solve!(p, solver)
            @test p.optval ≈ 19 atol=TOL
            @test evaluate(dotsort(x, [1, 2, 3, 4])) ≈ 19 atol=TOL
        end

        @testset "hinge loss atom" begin
            # TODO: @davidlizeng. We should finish this someday.
        end

        @testset "norm inf atom" begin
            x = Variable(3)
            p = minimize(norm_inf(x), [-2 <= x, x <= 1])
            @test vexity(p) == ConvexVexity()
            solve!(p, solver)
            @test p.optval ≈ 0 atol=TOL
            @test evaluate(norm_inf(x)) ≈ 0 atol=TOL
        end

        @testset "norm 1 atom" begin
            x = Variable(3)
            p = minimize(norm_1(x), [-2 <= x, x <= 1])
            @test vexity(p) == ConvexVexity()
            solve!(p, solver)
            @test p.optval ≈ 0 atol=TOL
            @test evaluate(norm_1(x)) ≈ 0 atol=TOL
        end

    end

    if can_solve_socp(solver)
        @testset "SOCP Atoms" begin

            @testset "norm 2 atom" begin
                x = Variable(2, 1)
                A = [1 2; 2 1; 3 4]
                b = [2; 3; 4]
                p = minimize(norm2(A * x + b))
                @test vexity(p) == ConvexVexity()
                solve!(p, solver)
                @test p.optval ≈ 0.64888 atol=TOL
                @test evaluate(norm2(A * x + b)) ≈ 0.64888 atol=TOL

                x = Variable(2, 1)
                A = [1 2; 2 1; 3 4]
                b = [2; 3; 4]
                lambda = 1
                p = minimize(norm2(A * x + b) + lambda * norm2(x), x >= 1)
                @test vexity(p) == ConvexVexity()
                solve!(p, solver)
                @test p.optval ≈ 14.9049 atol=TOL
                @test evaluate(norm2(A * x + b) + lambda * norm2(x)) ≈ 14.9049 atol=TOL

                x = Variable(2)

                p = minimize(norm2([x[1] + 2x[2] + 2; 2x[1] + x[2] + 3; 3x[1]+4x[2] + 4]) + lambda * norm2(x), x >= 1)
                @test vexity(p) == ConvexVexity()

                solve!(p, solver)
                @test p.optval ≈ 14.9049 atol=TOL
                @test evaluate(norm2(A * x + b) + lambda * norm2(x)) ≈ 14.9049 atol=TOL

                x = Variable(2, 1)
                A = [1 2; 2 1; 3 4]
                b = [2; 3; 4]
                lambda = 1
                p = minimize(norm2(A * x + b) + lambda * norm_1(x), x >= 1)
                @test vexity(p) == ConvexVexity()
                solve!(p, solver)
                @test p.optval ≈ 15.4907 atol=TOL
                @test evaluate(norm2(A * x + b) + lambda * norm_1(x)) ≈ 15.4907 atol=TOL
            end

            @testset "frobenius norm atom" begin
                m = Variable(4, 5)
                c = [m[3, 3] == 4, m >= 1]
                p = minimize(norm(vec(m), 2), c)
                @test vexity(p) == ConvexVexity()
                solve!(p, solver)
                @test p.optval ≈ sqrt(35) atol=TOL
                @test evaluate(norm(vec(m), 2)) ≈ sqrt(35) atol=TOL
            end

            @testset "quad over lin atom" begin
                x = Variable(3, 1)
                A = [2 -3 5; -2 9 -3; 5 -8 3]
                b = [-3; 9; 5]
                c = [3 2 4]
                d = -3
                p = minimize(quadoverlin(A*x + b, c*x + d))
                @test vexity(p) == ConvexVexity()
                solve!(p, solver)
                @test p.optval ≈ 17.7831 atol=TOL
                @test (evaluate(quadoverlin(A * x + b, c * x + d)))[1] ≈ 17.7831 atol=TOL
            end

            @testset "sum squares atom" begin
                x = Variable(2, 1)
                A = [1 2; 2 1; 3 4]
                b = [2; 3; 4]
                p = minimize(sumsquares(A*x + b))
                @test vexity(p) == ConvexVexity()
                solve!(p, solver)
                @test p.optval ≈ 0.42105 atol=TOL
                @test (evaluate(sumsquares(A * x + b)))[1] ≈ 0.42105 atol=TOL
            end

            @testset "square atom" begin
                x = Variable(2, 1)
                A = [1 2; 2 1; 3 4]
                b = [2; 3; 4]
                p = minimize(sum(square(A*x + b)))
                @test vexity(p) == ConvexVexity()
                solve!(p, solver)
                @test p.optval ≈ 0.42105 atol=TOL
                @test evaluate(sum(square(A * x + b))) ≈ 0.42105 atol=TOL

                x = Variable(2, 1)
                A = [1 2; 2 1; 3 4]
                b = [2; 3; 4]
                expr = A * x + b
                p = minimize(sum(dot(^)(expr,2))) # elementwise ^
                @test vexity(p) == ConvexVexity()
                solve!(p, solver)
                @test p.optval ≈ 0.42105 atol=TOL
                @test evaluate(sum(broadcast(^, expr, 2))) ≈ 0.42105 atol=TOL

                p = minimize(sum(dot(*)(expr, expr))) # elementwise *
                @test vexity(p) == ConvexVexity()
                solve!(p, solver)
                @test p.optval ≈ 0.42105 atol=TOL
                @test evaluate(sum((dot(*))(expr, expr))) ≈ 0.42105 atol=TOL
            end

            @testset "inv pos atom" begin
                x = Variable(4)
                p = minimize(sum(invpos(x)), invpos(x) < 2, x > 1, x == 2, 2 == x)
                @test vexity(p) == ConvexVexity()
                solve!(p, solver)
                @test p.optval ≈ 2 atol=TOL
                @test evaluate(sum(invpos(x))) ≈ 2 atol=TOL

                x = Variable(3)
                p = minimize(sum(dot(/)([3,6,9], x)), x<=3)
                solve!(p, solver)
                @test x.value ≈ fill(3.0, (3, 1)) atol=TOL
                @test p.optval ≈ 6 atol=TOL
                @test evaluate(sum((dot(/))([3, 6, 9], x))) ≈ 6 atol=TOL

                x = Variable()
                p = minimize(sum([3,6,9]/x), x<=3)
                solve!(p, solver)
                @test x.value ≈ 3 atol=TOL
                @test p.optval ≈ 6 atol=TOL
                @test evaluate(sum([3, 6, 9] / x)) ≈ 6 atol=TOL
            end

            @testset "geo mean atom" begin
                x = Variable(2)
                y = Variable(2)
                p = minimize(geomean(x, y), x >= 1, y >= 2)
                # not DCP compliant
                @test vexity(p) == ConcaveVexity()
                p = maximize(geomean(x, y), 1 < x, x < 2, y < 2)
                # Just gave it a vector as an objective, not okay
                @test_throws Exception solve!(p, solver)

                p = maximize(sum(geomean(x, y)), 1 < x, x < 2, y < 2)
                solve!(p, solver)
                @test p.optval ≈ 4 atol=TOL
                @test evaluate(sum(geomean(x, y))) ≈ 4 atol=TOL
            end

            @testset "sqrt atom" begin
                x = Variable()
                p = maximize(sqrt(x), 1 >= x)
            end

            @testset "quad form atom" begin
                x = Variable(3, 1)
                A = [0.8608 0.3131 0.5458; 0.3131 0.8584 0.5836; 0.5458 0.5836 1.5422]
                p = minimize(quadform(x, A), [x >= 1])
                @test vexity(p) == ConvexVexity()
                solve!(p, solver)
                @test p.optval ≈ 6.1464 atol=TOL
                @test (evaluate(quadform(x, A)))[1] ≈ 6.1464 atol=TOL

                x = Variable(3, 1)
                A = -1.0*[0.8608 0.3131 0.5458; 0.3131 0.8584 0.5836; 0.5458 0.5836 1.5422]
                c = [3 2 4]
                p = maximize(c*x , [quadform(x, A) >= -1])
                @test vexity(p) == ConvexVexity()
                solve!(p, solver)
                @test p.optval ≈ 3.7713 atol=TOL
                @test (evaluate(quadform(x, A)))[1] ≈ -1 atol=TOL
            end

            @testset "huber atom" begin
                x = Variable(3)
                p = minimize(sum(huber(x, 1)), x >= 2)
                @test vexity(p) == ConvexVexity()
                solve!(p, solver)
                @test p.optval ≈ 9 atol=TOL
                @test evaluate(sum(huber(x, 1))) ≈ 9 atol=TOL
            end

            @testset "rational norm atom" begin
                A = [1 2 3; -1 2 3]
                b = A * ones(3)
                x = Variable(3)
                p = minimize(norm(x, 4.5), [A * x == b])
                @test vexity(p) == ConvexVexity()
                # Solution is approximately x = [1, .93138, 1.04575]
                solve!(p, solver)
                @test p.optval ≈ 1.2717 atol=TOL
                @test evaluate(norm(x, 4.5)) ≈ 1.2717 atol=TOL
            end

            @testset "rational norm dual norm" begin
                v = [0.463339, 0.0216084, -2.07914, 0.99581, 0.889391]
                x = Variable(5)
                q = 1.379;  # q norm constraint that generates many inequalities
                qs = q / (q - 1);  # Conjugate to q
                p = minimize(x' * v)
                p.constraints += (norm(x, q) <= 1)
                @test vexity(p) == ConvexVexity()
                solve!(p, solver)  # Solution is -norm(v, q / (q - 1))
                @test p.optval ≈ -2.144087 atol=TOL
                @test sum(evaluate(x' * v)) ≈ -2.144087 atol=TOL
                @test evaluate(norm(x, q)) ≈ 1 atol=TOL
                @test sum(evaluate(x' * v)) ≈ -(sum(abs.(v) .^ qs) ^ (1 / qs)) atol=TOL
            end

            @testset "rational norm atom sum" begin
                A = [-0.719255  -0.229089
                    -1.33632   -1.37121
                    0.703447  -1.4482]
                b = [-1.82041, -1.67516, -0.866884]
                q = 1.5
                xvar = Variable(2)
                p = minimize(.5 * sumsquares(xvar) + norm(A * xvar - b, q))
                @test vexity(p) == ConvexVexity()
                solve!(p, solver)
                # Compute gradient, check it is zero(ish)
                x_opt = xvar.value
                margins = A * x_opt - b
                qs = q / (q - 1);  # Conjugate
                denom = sum(abs.(margins).^q)^(1/qs)
                g = x_opt + A' * (abs.(margins).^(q-1) .* sign.(margins)) / denom
                @test p.optval ≈ 1.7227 atol=TOL
                @test norm(g, 2) ^ 2 ≈ 0 atol=TOL
            end

            @testset "norm consistent with Base for matrix variables" begin
                A = randn(4, 4)
                x = Variable(4, 4)
                x.value = A
                # Matrix norm
                @test evaluate(opnorm(x)) ≈ opnorm(A) atol=TOL
                @test evaluate(opnorm(x, 1)) ≈ opnorm(A, 1) atol=TOL
                @test evaluate(opnorm(x, 2)) ≈ opnorm(A, 2) atol=TOL
                @test evaluate(opnorm(x, Inf)) ≈ opnorm(A, Inf) atol=TOL
                # Vector norm
                # TODO: Once the deprecation for norm on matrices is removed, remove the `vec` calls
                @test evaluate(norm(vec(x), 1)) ≈ norm(vec(A), 1) atol=TOL
                @test evaluate(norm(vec(x), 2)) ≈ norm(vec(A), 2) atol=TOL
                @test evaluate(norm(vec(x), 7)) ≈ norm(vec(A), 7) atol=TOL
                @test evaluate(norm(vec(x), Inf)) ≈ norm(vec(A), Inf) atol=TOL
            end
        end

        @testset "Fixed and freed variables" begin
            @testset "fix and free addition" begin
                x = Variable()
                y = Variable()

                p = minimize(x+y, x>=0, y>=0)
                solve!(p, solver)
                @test p.optval ≈ 0 atol=TOL

                y.value = 4
                fix!(y)
                solve!(p, solver)
                @test p.optval ≈ 4 atol=TOL

                free!(y)
                solve!(p, solver)
                @test p.optval ≈ 0 atol=TOL
            end

            @testset "fix multiplication" begin
                a = [1,2,3,2,1]
                x = Variable(length(a))
                gamma = Variable(Positive())
                fix!(gamma, 0.7)

                p = minimize(norm(x-a) + gamma*norm(x[1:end-1] - x[2:end]))
                solve!(p, solver)
                o1 = p.optval
                # x should be very close to a
                @test o1 ≈ 0.7 * norm(a[1:end - 1] - a[2:end]) atol=TOL
                # increase regularization
                fix!(gamma, 1.0)
                solve!(p, solver)
                o2 = p.optval
                # x should be very close to mean(a)
                @test o2 ≈ norm(a .- mean(a)) atol=TOL

                @test o1 <= o2
            end
        end
    end

    if can_solve_sdp(solver)
        # TODO: uncomment vexity checks once SDP on vars/constraints changes vexity of problem

        @testset "SDP Atoms" begin
            @testset "sdp variables" begin
                y = Variable((2,2), :Semidefinite)
                p = minimize(y[1,1])
                # @fact vexity(p) --> ConvexVexity()
                solve!(p, solver)
                @test p.optval ≈ 0 atol=TOL

                y = Variable((3,3), :Semidefinite)
                p = minimize(y[1,1], y[2,2]==1)
                # @fact vexity(p) --> ConvexVexity()
                solve!(p, solver)
                @test p.optval ≈ 0 atol=TOL

                # Solution is obtained as y[2,2] -> infinity
                # This test fails on Mosek. See
                # https://github.com/JuliaOpt/Mosek.jl/issues/29
                # y = Variable((2, 2), :Semidefinite)
                # p = minimize(y[1, 1], y[1, 2] == 1)
                # # @fact vexity(p) --> ConvexVexity()
                # solve!(p, solver)
                # @fact p.optval --> roughly(0, TOL)

                y = Semidefinite(3)
                p = minimize(sum(diag(y)), y[1, 1] == 1)
                # @fact vexity(p) --> ConvexVexity()
                solve!(p, solver)
                @test p.optval ≈ 1 atol=TOL

                y = Variable((3, 3), :Semidefinite)
                p = minimize(tr(y), y[2,1]<=4, y[2,2]>=3)
                # @fact vexity(p) --> ConvexVexity()
                solve!(p, solver)
                @test p.optval ≈ 3 atol=TOL

                x = Variable(Positive())
                y = Semidefinite(3)
                p = minimize(y[1, 2], y[2, 1] == 1)
                # @fact vexity(p) --> ConvexVexity()
                solve!(p, solver)
                @test p.optval ≈ 1 atol=TOL
            end

            @testset "sdp constraints" begin
                # This test fails on Mosek
                x = Variable(Positive())
                y = Variable((3, 3))
                p = minimize(x + y[1, 1], isposdef(y), x >= 1, y[2, 1] == 1)
                # @fact vexity(p) --> ConvexVexity()
                solve!(p, solver)
                @test p.optval ≈ 1 atol=TOL
            end

            @testset "nuclear norm atom" begin
                y = Semidefinite(3)
                p = minimize(nuclearnorm(y), y[2,1]<=4, y[2,2]>=3, y[3,3]<=2)
                @test vexity(p) == ConvexVexity()
                solve!(p, solver)
                @test p.optval ≈ 3 atol=TOL
                @test evaluate(nuclearnorm(y)) ≈ 3 atol=TOL
            end

            @testset "operator norm atom" begin
                y = Variable((3,3))
                p = minimize(opnorm(y), y[2,1]<=4, y[2,2]>=3, sum(y)>=12)
                @test vexity(p) == ConvexVexity()
                solve!(p, solver)
                @test p.optval ≈ 4 atol=TOL
                @test evaluate(opnorm(y)) ≈ 4 atol=TOL
            end

            @testset "sigma max atom" begin
                y = Variable((3,3))
                p = minimize(sigmamax(y), y[2,1]<=4, y[2,2]>=3, sum(y)>=12)
                @test vexity(p) == ConvexVexity()
                solve!(p, solver)
                @test p.optval ≈ 4 atol=TOL
                @test evaluate(sigmamax(y)) ≈ 4 atol=TOL
            end

            @testset "lambda max atom" begin
                y = Semidefinite(3)
                p = minimize(lambdamax(y), y[1,1]>=4)
                @test vexity(p) == ConvexVexity()
                solve!(p, solver)
                @test p.optval ≈ 4 atol=TOL
                @test evaluate(lambdamax(y)) ≈ 4 atol=TOL
            end

            @testset "lambda min atom" begin
                y = Semidefinite(3)
                p = maximize(lambdamin(y), tr(y)<=6)
                @test vexity(p) == ConvexVexity()
                solve!(p, solver)
                @test p.optval ≈ 2 atol=TOL
                @test evaluate(lambdamin(y)) ≈ 2 atol=TOL
            end

            @testset "matrix frac atom" begin
                x = [1, 2, 3]
                P = Variable(3, 3)
                p = minimize(matrixfrac(x, P), P <= 2*eye(3), P >= 0.5 * eye(3))
                @test vexity(p) == ConvexVexity()
                solve!(p, solver)
                @test p.optval ≈ 7 atol=TOL
                @test (evaluate(matrixfrac(x, P)))[1] ≈ 7 atol=TOL
            end

            @testset "matrix frac atom both arguments variable" begin
                x = Variable(3)
                P = Variable(3, 3)
                p = minimize(matrixfrac(x, P), lambdamax(P) <= 2, x[1] >= 1)
                @test vexity(p) == ConvexVexity()
                solve!(p, solver)
                @test p.optval ≈ 0.5 atol=TOL
                @test (evaluate(matrixfrac(x, P)))[1] ≈ 0.5 atol=TOL
            end

            @testset "sum largest eigs" begin
                x = Semidefinite(3)
                p = minimize(sumlargesteigs(x, 2), x >= 1)
                solve!(p, solver)
                @test p.optval ≈ 3 atol=TOL
                @test evaluate(x) ≈ ones(3, 3) atol=TOL

                x = Semidefinite(3)
                p = minimize(sumlargesteigs(x, 2), [x[i,:] >= i for i=1:3]...)
                solve!(p, solver)
                @test p.optval ≈ 8.4853 atol=TOL

                x1 = Semidefinite(3)
                p1 = minimize(lambdamax(x1), x1[1,1]>=4)
                solve!(p1, solver)

                x2 = Semidefinite(3)
                p2 = minimize(sumlargesteigs(x2, 1), x2[1,1]>=4)
                solve!(p2, solver)

                @test p1.optval ≈ p2.optval atol=TOL

                x1 = Semidefinite(3)
                p1 = minimize(lambdamax(x1), [x1[i,:] >= i for i=1:3]...)
                solve!(p1, solver)

                x2 = Semidefinite(3)
                p2 = minimize(sumlargesteigs(x2, 1), [x2[i,:] >= i for i=1:3]...)
                solve!(p2, solver)

                @test p1.optval ≈ p2.optval atol=TOL

                println(p1.optval)
            end

            @testset "kron atom" begin
                id = eye(4)
                X = Semidefinite(4)
                W = kron(id, X)
                p = maximize(tr(W), tr(X) ≤ 1)
                @test vexity(p) == AffineVexity()
                solve!(p, solver)
                @test p.optval ≈ 4 atol=TOL
            end

            @testset "Partial trace" begin
                A = Semidefinite(2)
                B = [1 0; 0 0]
                ρ = kron(B, A)
                constraints = [partialtrace(ρ, 1, [2; 2]) == [0.09942819 0.29923607; 0.29923607 0.90057181], ρ in :SDP]
                p = satisfy(constraints)
                solve!(p, solver)
                @test evaluate(ρ) ≈ [0.09942819 0.29923607 0 0; 0.299237 0.900572 0 0; 0 0 0 0; 0 0 0 0] atol=TOL
                @test evaluate(partialtrace(ρ, 1, [2; 2])) ≈ [1.0 0; 0 0] atol=TOL
            end
        end

        @testset "Optimization with complex variables" begin
            @testset "Real Variables with complex equality constraints" begin
                n = 10 # variable dimension (parameter)
                m = 5 # number of constraints (parameter)
                xo = rand(n)
                A = randn(m,n) + im*randn(m,n)
                b = A * xo
                x = Variable(n)
                p1 = minimize(sum(x), A*x == b, x>=0)
                solve!(p1, solver)
                x1 = x.value

                p2 = minimize(sum(x), real(A)*x == real(b), imag(A)*x==imag(b), x>=0)
                solve!(p2, solver)
                x2 = x.value
                @test x1 == x2
            end

            @testset "Complex Variable with complex equality constraints" begin
                n = 10 # variable dimension (parameter)
                m = 5 # number of constraints (parameter)
                xo = rand(n)+im*rand(n)
                A = randn(m,n) + im*randn(m,n)
                b = A * xo
                x = ComplexVariable(n)
                p1 = minimize(real(sum(x)), A*x == b, real(x)>=0, imag(x)>=0)
                solve!(p1, solver)
                x1 = x.value

                xr = Variable(n)
                xi = Variable(n)
                p2 = minimize(sum(xr), real(A)*xr-imag(A)*xi == real(b), imag(A)*xr+real(A)*xi == imag(b), xr>=0, xi>=0)
                solve!(p2, solver)
                #x2 = xr.value + im*xi.value
                real_diff = real(x1) - xr.value

                @test real_diff ≈ zeros(10, 1) atol=TOL
                imag_diff = imag(x1) - xi.value
                @test imag_diff ≈ zeros(10, 1) atol=TOL
                #@fact x1==x2 --> true
            end

            @testset "norm2 atom" begin
                a = 2+4im
                x = ComplexVariable()
                objective = norm2(a-x)
                c1 = real(x)>=0
                p = minimize(objective,c1)
                solve!(p, solver)
                @test p.optval ≈ 0 atol=TOL
                @test evaluate(objective) ≈ 0 atol=TOL
                real_diff = real(x.value) - real(a)
                imag_diff = imag(x.value) - imag(a)
                @test real_diff ≈ 0 atol=TOL
                @test imag_diff ≈ 0 atol=TOL
            end

            @testset "sumsquares atom" begin
                a = [2+4im;4+6im]
                x = ComplexVariable(2)
                objective = sumsquares(a-x)
                c1 = real(x)>=0
                p = minimize(objective,c1)
                solve!(p, solver)
                @test p.optval ≈ 0 atol=TOL
                @test evaluate(objective) ≈ zeros(1, 1) atol=TOL
                real_diff = real.(x.value) - real.(a)
                imag_diff = imag.(x.value) - imag.(a)
                @test real_diff ≈ zeros(2, 1) atol=TOL
                @test imag_diff ≈ zeros(2, 1) atol=TOL
            end

            @testset "abs atom" begin
                a = [5-4im]
                x = ComplexVariable()
                objective = abs(a-x)
                c1 = real(x)>=0
                p = minimize(objective,c1)
                solve!(p, solver)
                @test p.optval ≈ 0 atol=TOL
                @test evaluate(objective) ≈ zeros(1) atol=TOL
                real_diff = real(x.value) .- real(a)
                imag_diff = imag(x.value) .- imag(a)
                @test real_diff ≈ zeros(1) atol=TOL
                @test imag_diff ≈ zeros(1) atol=TOL
            end

            @testset "Complex Semidefinite constraint" begin
                n = 10
                A = rand(n,n) + im*rand(n,n)
                A = A + A' # now A is hermitian
                x = ComplexVariable(n,n)
                objective = sumsquares(A - x)
                c1 = x in :SDP
                p = minimize(objective, c1)
                solve!(p, solver)
                # test that X is approximately equal to posA:
                l,v = eigen(A)
                posA = v*Diagonal(max.(l,0))*v'

                real_diff = real.(x.value) - real.(posA)
                imag_diff = imag.(x.value) - imag.(posA)
                @test real_diff ≈ zeros(n, n) atol=TOL
                @test imag_diff ≈ zeros(n, n) atol=TOL
            end
        end
    end

    if can_solve_exp(solver)
        @testset "Exp Atoms" begin
            @testset "exp atom" begin
                y = Variable()
                p = minimize(exp(y), y>=0)
                @test vexity(p) == ConvexVexity()
                solve!(p, solver)
                @test p.optval ≈ 1 atol=TOL
                @test evaluate(exp(y)) ≈ 1 atol=TOL

                y = Variable()
                p = minimize(exp(y), y>=1)
                @test vexity(p) == ConvexVexity()
                solve!(p, solver)
                @test p.optval ≈ exp(1) atol=TOL
                @test evaluate(exp(y)) ≈ exp(1) atol=TOL

                y = Variable(5)
                p = minimize(sum(exp(y)), y>=0)
                @test vexity(p) == ConvexVexity()
                solve!(p, solver)
                @test p.optval ≈ 5 atol=TOL
                @test evaluate(sum(exp(y))) ≈ 5 atol=TOL

                y = Variable(5)
                p = minimize(sum(exp(y)), y>=0)
                @test vexity(p) == ConvexVexity()
                solve!(p, solver)
                @test p.optval ≈ 5 atol=TOL
            end

            @testset "log atom" begin
                y = Variable()
                p = maximize(log(y), y<=1)
                @test vexity(p) == ConvexVexity()
                solve!(p, solver)
                @test p.optval ≈ 0 atol=TOL

                y = Variable()
                p = maximize(log(y), y<=2)
                @test vexity(p) == ConvexVexity()
                solve!(p, solver)
                @test p.optval ≈ log(2) atol=TOL

                y = Variable()
                p = maximize(log(y), [y<=2, exp(y)<=10])
                @test vexity(p) == ConvexVexity()
                solve!(p, solver)
                @test p.optval ≈ log(2) atol=TOL
            end

            @testset "log sum exp atom" begin
                y = Variable(5)
                p = minimize(logsumexp(y), y>=1)
                @test vexity(p) == ConvexVexity()
                solve!(p, solver)
                @test p.optval ≈ log(exp(1) * 5) atol=TOL
            end

            @testset "logistic loss atom" begin
                y = Variable(5)
                p = minimize(logisticloss(y), y>=1)
                @test vexity(p) == ConvexVexity()
                solve!(p, solver)
                @test p.optval ≈ log(exp(1) + 1) * 5 atol=TOL
            end

            @testset "entropy atom" begin
                y = Variable(5, Positive())
                p = maximize(entropy(y), sum(y)<=1)
                @test vexity(p) == ConvexVexity()
                solve!(p, solver)
                @test p.optval ≈ -(log(1 / 5)) atol=TOL
            end

            @testset "relative entropy atom" begin
                x = Variable(1)
                y = Variable(1)
                # x log (x/y)
                p = minimize(relative_entropy(x,y), y==1, x >= 2)
                @test vexity(p) == ConvexVexity()
                solve!(p, solver)
                @test p.optval ≈ 2 * log(2) atol=TOL
            end

            @testset "log perspective atom" begin
                x = Variable(1)
                y = Variable(1)
                # y log (x/y)
                p = maximize(log_perspective(x,y), y==5, x <= 10)
                @test vexity(p) == ConvexVexity()
                solve!(p, solver)
                @test p.optval ≈ 5 * log(2) atol=TOL
            end

        end

    end

    if can_solve_sdp(solver) && can_solve_exp(solver)
        @testset "SDP and Exp Atoms" begin
            tol = 1e-2
            @testset "log det atom" begin
                x = Variable(2, 2)
                p = maximize(logdet(x), [x[1, 1] == 1, x[2, 2] == 1])
                @test vexity(p) == ConvexVexity()
                solve!(p, solver)
                @test p.optval ≈ 0 atol=tol
                @test evaluate(logdet(x)) ≈ 0 atol=tol
            end
        end
    end

    if can_solve_mip(solver)
        mip_solver = !can_solve_mip(solver) ? GLPKSolverMIP() : solver

        @testset "Mixed Integer Programs" begin

            @testset "lp fallback interface" begin
                x = Variable()
                p = minimize(x, x>=4.3)
                @test vexity(p) == AffineVexity()
                solve!(p, solver)
                @test p.optval ≈ 4.3 atol=TOL

                x = Variable(2)
                p = minimize(norm(x,1), x[1]>=4.3)
                @test vexity(p) == ConvexVexity()
                solve!(p, solver)
                @test p.optval ≈ 4.3 atol=TOL
            end

            @testset "integer variables" begin
                x = Variable(:Int)
                p = minimize(x, x>=4.3)
                @test vexity(p) == AffineVexity()
                solve!(p, mip_solver)
                @test p.optval ≈ 5 atol=TOL

                x = Variable(2, :Int)
                p = minimize(sum(x), x>=4.3)
                @test vexity(p) == AffineVexity()
                solve!(p, mip_solver)
                @test p.optval ≈ 10 atol=TOL

                x = Variable(:Int)
                y = Variable()
                p = minimize(sum(x + y), x>=4.3, y>=7)
                @test vexity(p) == AffineVexity()
                solve!(p, mip_solver)
                @test p.optval ≈ 12 atol=TOL

                x = Variable(2, :Int)
                p = minimize(norm(x, 1), x[1]>=4.3)
                @test vexity(p) == ConvexVexity()
                solve!(p, mip_solver)
                @test p.optval ≈ 5 atol=TOL

                x = Variable(2, :Int)
                p = minimize(sum(x), x[1]>=4.3, x>=0)
                @test vexity(p) == AffineVexity()
                solve!(p, mip_solver)
                @test p.optval ≈ 5 atol=TOL

                x = Variable(2, :Int)
                p = minimize(sum(x), x>=.5)
                @test vexity(p) == AffineVexity()
                solve!(p, mip_solver)
                @test p.optval ≈ 2 atol=TOL
            end

            @testset "binary variables" begin
                x = Variable(2, :Bin)
                p = minimize(sum(x), x>=.5)
                @test vexity(p) == AffineVexity()
                solve!(p, mip_solver)
                @test p.optval ≈ 2 atol=TOL

                x = Variable(2, :Bin)
                p = minimize(sum(x), x[1]>=.5, x>=0)
                @test vexity(p) == AffineVexity()
                solve!(p, mip_solver)
                @test p.optval ≈ 1 atol=TOL
            end
        end
    end
end

# end
# tests = ["test_utilities.jl",
#          "test_affine.jl",
#          "test_lp.jl"]
# tests_socp = ["test_socp.jl","test_params.jl"]
# tests_sdp = ["test_sdp.jl"]
# tests_exp = ["test_exp.jl"]
# tests_int = ["test_int.jl"]
# tests_exp_and_sdp = ["test_exp_and_sdp.jl"]
# tests_complex = ["test_complex.jl"]

# println("Running tests:")

# for curtest in tests
#     @info " Test: $(curtest)"
#     include(curtest)
# end

# if can_solve_socp(solver)
#     for curtest in tests_socp
#         @info " Test: $(curtest)"
#         include(curtest)
#     end
# end

# if can_solve_sdp(solver)
#     for curtest in tests_sdp
#         @info " Test: $(curtest)"
#         include(curtest)
#     end
# end

# if can_solve_exp(solver)
#     for curtest in tests_exp
#         @info " Test: $(curtest)"
#         include(curtest)
#     end
# end

# if can_solve_sdp(solver) && can_solve_exp(solver)
#     for curtest in tests_exp_and_sdp
#         @info " Test: $(curtest)"
#         include(curtest)
#     end
# end

# if can_solve_mip(solver)
#     for curtest in tests_int
#         @info " Test: $(curtest)"
#         include(curtest)
#     end
# end

# if can_solve_sdp(solver)
#     for curtest in tests_complex
#         @info " Test: $(curtest)"
#         include(curtest)
#     end
# end

# end
