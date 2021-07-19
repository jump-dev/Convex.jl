# TODO: uncomment vexity checks once SDP on vars/constraints changes vexity of problem
@add_problem sdp function sdp_sdp_variables(handle_problem!, ::Val{test}, atol, rtol, ::Type{T}) where {T, test}
    y = Variable((2,2), :Semidefinite)
    p = minimize(y[1,1]; numeric_type = T)

    # @fact vexity(p) --> ConvexVexity()
    handle_problem!(p)
    if test
        @test p.optval ≈ 0 atol=atol rtol=rtol
    end

    y = Variable((3,3), :Semidefinite)
    p = minimize(y[1,1], y[2,2]==1; numeric_type = T)

    # @fact vexity(p) --> ConvexVexity()
    handle_problem!(p)
    if test
        @test p.optval ≈ 0 atol=atol rtol=rtol
    end

    # Solution is obtained as y[2,2] -> infinity
    # This test fails on Mosek. See
    # https://github.com/JuliaOpt/Mosek.jl/issues/29
    # y = Variable((2, 2), :Semidefinite)
    # p = minimize(y[1, 1], y[1, 2] == 1; numeric_type = T)

    # # @fact vexity(p) --> ConvexVexity()
    # handle_problem!(p)
    # @fact p.optval --> roughly(0, atol)

    y = Semidefinite(3)
    p = minimize(sum(diag(y)), y[1, 1] == 1; numeric_type = T)

    # @fact vexity(p) --> ConvexVexity()
    handle_problem!(p)
    if test
        @test p.optval ≈ 1 atol=atol rtol=rtol
    end

    y = Variable((3, 3), :Semidefinite)
    p = minimize(tr(y), y[2,1]<=4, y[2,2]>=3; numeric_type = T)

    # @fact vexity(p) --> ConvexVexity()
    handle_problem!(p)
    if test
        @test p.optval ≈ 3 atol=atol rtol=rtol
    end

    x = Variable(Positive())
    y = Semidefinite(3)
    p = minimize(y[1, 2], y[2, 1] == 1; numeric_type = T)

    # @fact vexity(p) --> ConvexVexity()
    handle_problem!(p)
    if test
        @test p.optval ≈ 1 atol=atol rtol=rtol
    end
end

@add_problem sdp function sdp_sdp_constraints(handle_problem!, ::Val{test}, atol, rtol, ::Type{T}) where {T, test}
    x = Variable(Positive())
    y = Variable((3, 3))
    p = minimize(x + y[1, 1], isposdef(y), x >= 1, y[2, 1] == 1; numeric_type = T)

    # @fact vexity(p) --> ConvexVexity()
    handle_problem!(p)
    if test
        @test p.optval ≈ 1 atol=atol rtol=rtol
    end
end

@add_problem sdp function sdp_nuclear_norm_atom(handle_problem!, ::Val{test}, atol, rtol, ::Type{T}) where {T, test}
    y = Semidefinite(3)
    p = minimize(nuclearnorm(y), y[2,1]<=4, y[2,2]>=3, y[3,3]<=2; numeric_type = T)

    if test
        @test vexity(p) == ConvexVexity()
    end
    handle_problem!(p)
    if test
        @test p.optval ≈ 3 atol=atol rtol=rtol
        @test evaluate(nuclearnorm(y)) ≈ 3 atol=atol rtol=rtol
    end
end

@add_problem sdp function sdp_nuclear_norm_atom_complex(handle_problem!, ::Val{test}, atol, rtol, ::Type{T}) where {T, test}
    A = [1 2im 3 4; 4im 3im 2 1; 4 5 6 7]
    y = ComplexVariable(3, 4)
    p = minimize(nuclearnorm(y), y == A; numeric_type = T)

    if test
        @test vexity(p) == ConvexVexity()
    end
    handle_problem!(p)
    if test
        @test p.optval ≈ sum(svdvals(A)) atol=atol rtol=rtol
        @test evaluate(nuclearnorm(y)) ≈ sum(svdvals(A)) atol=atol rtol=rtol
    end
end

@add_problem sdp function sdp_operator_norm_atom(handle_problem!, ::Val{test}, atol, rtol, ::Type{T}) where {T, test}
    y = Variable((3,3))
    p = minimize(opnorm(y), y[2,1]<=4, y[2,2]>=3, sum(y)>=12; numeric_type = T)

    if test
        @test vexity(p) == ConvexVexity()
    end
    handle_problem!(p)
    if test
        @test p.optval ≈ 4 atol=atol rtol=rtol
        @test evaluate(opnorm(y)) ≈ 4 atol=atol rtol=rtol
    end
end

@add_problem sdp function sdp_operator_norm_atom_complex(handle_problem!, ::Val{test}, atol, rtol, ::Type{T}) where {T, test}
    A = [1 2im 3 4; 4im 3im 2 1; 4 5 6 7]
    y = ComplexVariable(3, 4)
    p = minimize(opnorm(y), y == A; numeric_type = T)

    if test
        @test vexity(p) == ConvexVexity()
    end
    handle_problem!(p)
    if test
        @test p.optval ≈ maximum(svdvals(A)) atol=atol rtol=rtol
        @test evaluate(opnorm(y)) ≈ maximum(svdvals(A)) atol=atol rtol=rtol
    end
end

@add_problem sdp function sdp_sigma_max_atom(handle_problem!, ::Val{test}, atol, rtol, ::Type{T}) where {T, test}
    y = Variable((3,3))
    p = minimize(sigmamax(y), y[2,1]<=4, y[2,2]>=3, sum(y)>=12; numeric_type = T)

    if test
        @test vexity(p) == ConvexVexity()
    end
    handle_problem!(p)
    if test
        @test p.optval ≈ 4 atol=atol rtol=rtol
        @test evaluate(sigmamax(y)) ≈ 4 atol=atol rtol=rtol
    end
end

@add_problem sdp function sdp_dual_lambda_max_atom(handle_problem!, ::Val{test}, atol, rtol, ::Type{T}) where {T, test}
    y = Semidefinite(3)
    p = minimize(eigmax(y), y[1,1]>=4; numeric_type = T)

    if test
        @test vexity(p) == ConvexVexity()
    end
    handle_problem!(p)
    if test
        @test p.optval ≈ 4 atol=atol rtol=rtol
        @test evaluate(eigmax(y)) ≈ 4 atol=atol rtol=rtol
    end

    # https://github.com/jump-dev/Convex.jl/issues/337
    x = ComplexVariable(2, 2)
    p = minimize( eigmax(x), [ x[1,2] == im, x[2,2] == 1.0, x ⪰ - eye(2) ]; numeric_type = T)
    handle_problem!(p)
    if test
        @test p.optval ≈ 1.5 atol=atol rtol=rtol
        @test p.constraints[1].dual ≈ im atol=atol rtol=rtol
        @test p.constraints[2].dual ≈ 0.75 atol=atol rtol=rtol
    end
end

@add_problem sdp function sdp_lambda_min_atom(handle_problem!, ::Val{test}, atol, rtol, ::Type{T}) where {T, test}
    y = Semidefinite(3)
    p = maximize(eigmin(y), tr(y)<=6; numeric_type = T)

    if test
        @test vexity(p) == ConvexVexity()
    end
    handle_problem!(p)
    if test
        @test p.optval ≈ 2 atol=atol rtol=rtol
        @test evaluate(eigmin(y)) ≈ 2 atol=atol rtol=rtol
    end
end

@add_problem sdp function sdp_matrix_frac_atom(handle_problem!, ::Val{test}, atol, rtol, ::Type{T}) where {T, test}
    x = [1, 2, 3]
    P = Variable(3, 3)
    p = minimize(matrixfrac(x, P), P <= 2*eye(3), P >= 0.5 * eye(3); numeric_type = T)

    if test
        @test vexity(p) == ConvexVexity()
    end
    handle_problem!(p)
    if test
        @test p.optval ≈ 7 atol=atol rtol=rtol
        @test (evaluate(matrixfrac(x, P)))[1] ≈ 7 atol=atol rtol=rtol
    end
end

@add_problem sdp function sdp_matrix_frac_atom_both_arguments_variable(handle_problem!, ::Val{test}, atol, rtol, ::Type{T}) where {T, test}
    x = Variable(3)
    P = Variable(3, 3)
    p = minimize(matrixfrac(x, P), eigmax(P) <= 2, x[1] >= 1; numeric_type = T)

    if test
        @test vexity(p) == ConvexVexity()
    end
    handle_problem!(p)
    if test
        @test p.optval ≈ 0.5 atol=atol rtol=rtol
        @test (evaluate(matrixfrac(x, P)))[1] ≈ 0.5 atol=atol rtol=rtol
    end
end

@add_problem sdp function sdp_sum_largest_eigs(handle_problem!, ::Val{test}, atol, rtol, ::Type{T}) where {T, test}
    x = Semidefinite(3)
    p = minimize(sumlargesteigs(x, 2), x >= 1; numeric_type = T)

    handle_problem!(p)

    if test
        @test p.optval ≈ 3 atol=atol rtol=rtol
        @test evaluate(x) ≈ ones(3, 3) atol=atol rtol=rtol
    end

    x = Semidefinite(3)
    p = minimize(sumlargesteigs(x, 2) + sumlargesteigs(x,0), [x[i,:] >= i for i=1:3]...; numeric_type = T)

    handle_problem!(p)

    if test
        @test p.optval ≈ 8.4853 atol=atol rtol=rtol
    end

    A = [1   -2im  3   4
         2im  7    3im 5
         3   -3im  2   9
         4    5    9   4]

    x = ComplexVariable(4, 4)
    p = minimize(sumlargesteigs(x, 3), x == A; numeric_type = T)

    handle_problem!(p)

    if test
        @test p.optval ≈ sum(eigvals(A)[2:end]) atol=atol rtol=rtol
    end

    x1 = Semidefinite(3)
    p1 = minimize(eigmax(x1), x1[1,1]>=4; numeric_type = T)

    handle_problem!(p1)

    x2 = Semidefinite(3)
    p2 = minimize(sumlargesteigs(x2, 1), x2[1,1]>=4; numeric_type = T)

    handle_problem!(p2)

    if test
        @test p1.optval ≈ p2.optval atol=atol rtol=rtol
    end

    x1 = Semidefinite(3)
    p1 = minimize(eigmax(x1), [x1[i,:] >= i for i=1:3]...; numeric_type = T)

    handle_problem!(p1)

    x2 = Semidefinite(3)
    p2 = minimize(sumlargesteigs(x2, 1), [x2[i,:] >= i for i=1:3]...; numeric_type = T)

    handle_problem!(p2)

    if test
        @test p1.optval ≈ p2.optval atol=atol rtol=rtol
    end

end

@add_problem sdp function sdp_kron_atom(handle_problem!, ::Val{test}, atol, rtol, ::Type{T}) where {T, test}
    id = eye(4)
    X = Semidefinite(4)
    W = kron(id, X)
    p = maximize(tr(W), tr(X) ≤ 1; numeric_type = T)

    if test
        @test vexity(p) == AffineVexity()
    end
    handle_problem!(p)
    if test
        @test p.optval ≈ 4 atol=atol rtol=rtol
    end
end

@add_problem sdp function sdp_Partial_trace(handle_problem!, ::Val{test}, atol, rtol, ::Type{T}) where {T, test}
    A = Semidefinite(2)
    B = [1 0; 0 0]
    ρ = kron(B, A)
    constraints = [partialtrace(ρ, 1, [2; 2]) == [0.09942819 0.29923607; 0.29923607 0.90057181], ρ in :SDP]
    p = satisfy(constraints; numeric_type = T)

    handle_problem!(p)
    if test
        @test evaluate(ρ) ≈ [0.09942819 0.29923607 0 0; 0.299237 0.900572 0 0; 0 0 0 0; 0 0 0 0] atol=atol rtol=rtol
        @test evaluate(partialtrace(ρ, 1, [2; 2])) ≈ [0.09942819 0.29923607; 0.29923607 0.90057181] atol=atol rtol=rtol
    end

    function rand_normalized(n)
        A = 5*randn(n, n) + im*5*randn(n, n)
        A / tr(A)
    end

    As = [ rand_normalized(3) for _ = 1:5]
    Bs = [ rand_normalized(2) for _ = 1:5]
    p = rand(5)

    AB = sum(i -> p[i]*kron(As[i],Bs[i]), 1:5)
    if test
        @test partialtrace(AB, 2, [3, 2]) ≈ sum( p .* As ) atol=atol rtol=rtol
        @test partialtrace(AB, 1, [3, 2]) ≈ sum( p .* Bs ) atol=atol rtol=rtol
    end

    A, B, C = rand(5,5), rand(4,4), rand(3,3)
    ABC = kron(kron(A, B), C)
    if test
        @test kron(A,B)*tr(C) ≈ partialtrace(ABC, 3, [5, 4, 3]) atol=atol rtol=rtol
    end

    # Test 281
    A = rand(6,6)
    expr = partialtrace(Constant(A), 1, [2, 3])
    if test
        @test size(expr) == size(evaluate(expr))

        @test_throws ArgumentError partialtrace(rand(6, 6), 3, [2, 3])
        @test_throws ArgumentError partialtrace(rand(6, 6), 1, [2, 4])
        @test_throws ArgumentError partialtrace(rand(3, 4), 1, [2, 3])
    end
end


@add_problem sdp function sdp_Real_Variables_with_complex_equality_constraints(handle_problem!, ::Val{test}, atol, rtol, ::Type{T}) where {T, test}
    n = 10 # variable dimension (parameter)
    m = 5 # number of constraints (parameter)
    xo = rand(n)
    A = randn(m,n) + im*randn(m,n)
    b = A * xo
    x = Variable(n)
    p1 = minimize(sum(x), A*x == b, x>=0; numeric_type = T)

    handle_problem!(p1)

    if test
        x1 = copy(evaluate(x))
    end

    p2 = minimize(sum(x), real(A)*x == real(b), imag(A)*x==imag(b), x>=0; numeric_type = T)

    handle_problem!(p2)
    if test
        x2 = evaluate(x)
        @test x1 == x2
    end
end

@add_problem sdp function sdp_Complex_Variable_with_complex_equality_constraints(handle_problem!, ::Val{test}, atol, rtol, ::Type{T}) where {T, test}
    n = 10 # variable dimension (parameter)
    m = 5 # number of constraints (parameter)
    xo = rand(n)+im*rand(n)
    A = randn(m,n) + im*randn(m,n)
    b = A * xo
    x = ComplexVariable(n)
    p1 = minimize(real(sum(x)), A*x == b, real(x)>=0, imag(x)>=0; numeric_type = T)

    handle_problem!(p1)

    if test
        x1 = evaluate(x)
    end

    xr = Variable(n)
    xi = Variable(n)
    p2 = minimize(sum(xr), real(A)*xr-imag(A)*xi == real(b), imag(A)*xr+real(A)*xi == imag(b), xr>=0, xi>=0; numeric_type = T)

    handle_problem!(p2)

    if test
        real_diff = real(x1) - evaluate(xr)
        @test real_diff ≈ zeros(10, 1) atol=atol rtol=rtol

        imag_diff = imag(x1) - evaluate(xi)
        @test imag_diff ≈ zeros(10, 1) atol=atol rtol=rtol
    end
end

@add_problem sdp function sdp_Issue_198(handle_problem!, ::Val{test}, atol, rtol, ::Type{T}) where {T, test}
    ρ = HermitianSemidefinite(2)
    constraints = [ρ == [ 1. 0.; 0.  1.]]
    p = satisfy(constraints; numeric_type = T)

    handle_problem!(p)
    if test
        @test p.status == MOI.OPTIMAL
        @test evaluate(ρ) ≈ [ 1. 0.; 0.  1.] atol=atol rtol=rtol
        @test p.optval ≈ 0 atol=atol rtol=rtol
    end
end

@add_problem sdp function sdp_socp_norm2_atom(handle_problem!, ::Val{test}, atol, rtol, ::Type{T}) where {T, test}
    a = 2+4im
    x = ComplexVariable()
    objective = norm2(a-x)
    c1 = real(x)>=0
    p = minimize(objective,c1; numeric_type = T)

    handle_problem!(p)
    if test
        @test p.optval ≈ 0 atol=atol rtol=rtol
        @test evaluate(objective) ≈ 0 atol=atol rtol=rtol

        real_diff = real(evaluate(x)) - real(a)
        imag_diff = imag(evaluate(x)) - imag(a)
        @test real_diff ≈ 0 atol=atol rtol=rtol
        @test imag_diff ≈ 0 atol=atol rtol=rtol
    end
end

@add_problem sdp function sdp_socp_sumsquares_atom(handle_problem!, ::Val{test}, atol, rtol, ::Type{T}) where {T, test}
    a = [2+4im;4+6im]
    x = ComplexVariable(2)
    objective = sumsquares(a-x)
    c1 = real(x)>=0
    p = minimize(objective,c1; numeric_type = T)

    handle_problem!(p)
    if test
        @test p.optval ≈ 0 atol=atol rtol=rtol
        @test evaluate(objective) ≈ 0.0 atol=atol rtol=rtol

        real_diff = real.(evaluate(x)) - real.(a)
        imag_diff = imag.(evaluate(x)) - imag.(a)
        @test real_diff ≈ zeros(2, 1) atol=atol rtol=rtol
        @test imag_diff ≈ zeros(2, 1) atol=atol rtol=rtol
    end
end

@add_problem sdp function sdp_socp_abs_atom(handle_problem!, ::Val{test}, atol, rtol, ::Type{T}) where {T, test}
    a = 5-4im
    x = ComplexVariable()
    objective = abs(a-x)
    c1 = real(x)>=0
    p = minimize(objective,c1; numeric_type = T)

    handle_problem!(p)
    if test
        @test p.optval ≈ 0 atol=atol rtol=rtol
        @test evaluate(objective) ≈ 0.0 atol=atol rtol=rtol

        real_diff = real(evaluate(x)) .- real(a)
        imag_diff = imag(evaluate(x)) .- imag(a)
        @test real_diff ≈ 0.0 atol=atol rtol=rtol
        @test imag_diff ≈ 0.0 atol=atol rtol=rtol
    end
end

@add_problem sdp function sdp_Complex_Semidefinite_constraint(handle_problem!, ::Val{test}, atol, rtol, ::Type{T}) where {T, test}
    n = 10
    A = rand(n,n) + im*rand(n,n)
    A = A + A' # now A is hermitian
    x = ComplexVariable(n,n)
    objective = sumsquares(A - x)
    c1 = x in :SDP
    p = minimize(objective, c1; numeric_type = T)

    handle_problem!(p)
    # test that X is approximately equal to posA:
    l,v = eigen(A)
    posA = v*Diagonal(max.(l,0))*v'


    if test
        real_diff = real.(evaluate(x)) - real.(posA)
        imag_diff = imag.(evaluate(x)) - imag.(posA)
        @test real_diff ≈ zeros(n, n) atol=atol rtol=rtol
        @test imag_diff ≈ zeros(n, n) atol=atol rtol=rtol
    end
end

@add_problem sdp function sdp_geom_mean_hypocone_argcheck(handle_problem!, ::Val{test}, atol, rtol, ::Type{T}) where {T, test}
    if test
        @test_throws DimensionMismatch GeomMeanHypoCone(   zeros(2,3),    zeros(2,3), 1//2)
        @test_throws DimensionMismatch GeomMeanHypoCone(   zeros(2,3), Variable(2,3), 1//2)
        @test_throws DimensionMismatch GeomMeanHypoCone(Variable(2,3),    zeros(2,3), 1//2)
        @test_throws DimensionMismatch GeomMeanHypoCone(Variable(2,3), Variable(2,3), 1//2)

        @test_throws DimensionMismatch GeomMeanHypoCone(   zeros(2,2),    zeros(3,3), 1//2)
        @test_throws DimensionMismatch GeomMeanHypoCone(   zeros(2,2), Variable(3,3), 1//2)
        @test_throws DimensionMismatch GeomMeanHypoCone(Variable(2,2),    zeros(3,3), 1//2)
        @test_throws DimensionMismatch GeomMeanHypoCone(Variable(2,2), Variable(2,3), 1//2)

        @test_throws DimensionMismatch Variable(2,2) in GeomMeanHypoCone(Variable(3,3), Variable(3,3), 1//2)

        @test_throws DomainError GeomMeanHypoCone(Variable(3,3), Variable(3,3), -1//2)
        @test_throws DomainError GeomMeanHypoCone(Variable(3,3), Variable(3,3), 3//2)
    end
end

@add_problem sdp function sdp_geom_mean_epicone_argcheck(handle_problem!, ::Val{test}, atol, rtol, ::Type{T}) where {T, test}
    if test
        @test_throws DimensionMismatch GeomMeanEpiCone(   zeros(2,3),    zeros(2,3), -1//2)
        @test_throws DimensionMismatch GeomMeanEpiCone(   zeros(2,3), Variable(2,3), -1//2)
        @test_throws DimensionMismatch GeomMeanEpiCone(Variable(2,3),    zeros(2,3), -1//2)
        @test_throws DimensionMismatch GeomMeanEpiCone(Variable(2,3), Variable(2,3), -1//2)

        @test_throws DimensionMismatch GeomMeanEpiCone(   zeros(2,2),    zeros(3,3), -1//2)
        @test_throws DimensionMismatch GeomMeanEpiCone(   zeros(2,2), Variable(3,3), -1//2)
        @test_throws DimensionMismatch GeomMeanEpiCone(Variable(2,2),    zeros(3,3), -1//2)
        @test_throws DimensionMismatch GeomMeanEpiCone(Variable(2,2), Variable(2,3), -1//2)

        @test_throws DimensionMismatch Variable(2,2) in GeomMeanEpiCone(Variable(3,3), Variable(3,3), -1//2)
        @test_throws DomainError GeomMeanEpiCone(Variable(3,3), Variable(3,3), -3//2)
        @test_throws DomainError GeomMeanEpiCone(Variable(3,3), Variable(3,3), 1//2)
        @test_throws DomainError GeomMeanEpiCone(Variable(3,3), Variable(3,3), 5//2)
    end
end

@add_problem sdp function sdp_relative_entropy_argcheck(handle_problem!, ::Val{test}, atol, rtol, ::Type{T}) where {T, test}
    if test
        @test_throws DimensionMismatch RelativeEntropyEpiCone(   zeros(2,3),    zeros(2,3), 3, 3, eye(2))
        @test_throws DimensionMismatch RelativeEntropyEpiCone(   zeros(2,3), Variable(2,3), 3, 3, eye(2))
        @test_throws DimensionMismatch RelativeEntropyEpiCone(Variable(2,3),    zeros(2,3), 3, 3, eye(2))
        @test_throws DimensionMismatch RelativeEntropyEpiCone(Variable(2,3), Variable(2,3), 3, 3, eye(2))

        @test_throws DimensionMismatch RelativeEntropyEpiCone(   zeros(2,2),    zeros(3,3), 3, 3)
        @test_throws DimensionMismatch RelativeEntropyEpiCone(   zeros(2,2), Variable(3,3), 3, 3)
        @test_throws DimensionMismatch RelativeEntropyEpiCone(Variable(2,2),    zeros(3,3), 3, 3)
        @test_throws DimensionMismatch RelativeEntropyEpiCone(Variable(2,2), Variable(2,3), 3, 3)

        @test_throws DimensionMismatch RelativeEntropyEpiCone(   zeros(3,3),    zeros(3,3), 3, 3, zeros(2, 2))
        @test_throws DimensionMismatch RelativeEntropyEpiCone(   zeros(3,3), Variable(3,3), 3, 3, zeros(2, 2))
        @test_throws DimensionMismatch RelativeEntropyEpiCone(Variable(3,3),    zeros(3,3), 3, 3, zeros(2, 2))
        @test_throws DimensionMismatch RelativeEntropyEpiCone(Variable(3,3), Variable(2,3), 3, 3, zeros(2, 2))

        @test_throws DimensionMismatch RelativeEntropyEpiCone(   zeros(3,3),    zeros(3,3), 3, 3, zeros(2, 3))
        @test_throws DimensionMismatch RelativeEntropyEpiCone(   zeros(3,3), Variable(3,3), 3, 3, zeros(2, 3))
        @test_throws DimensionMismatch RelativeEntropyEpiCone(Variable(3,3),    zeros(3,3), 3, 3, zeros(2, 3))
        @test_throws DimensionMismatch RelativeEntropyEpiCone(Variable(3,3), Variable(2,3), 3, 3, zeros(2, 3))

        @test_throws DimensionMismatch RelativeEntropyEpiCone(   zeros(3,3),    zeros(3,3), 3, 3, zeros(2))
        @test_throws DimensionMismatch RelativeEntropyEpiCone(   zeros(3,3), Variable(3,3), 3, 3, zeros(2))
        @test_throws DimensionMismatch RelativeEntropyEpiCone(Variable(3,3),    zeros(3,3), 3, 3, zeros(2))
        @test_throws DimensionMismatch RelativeEntropyEpiCone(Variable(3,3), Variable(2,3), 3, 3, zeros(2))

        @test_throws DimensionMismatch Variable(2,2) in RelativeEntropyEpiCone(Variable(3,3), Variable(3,3), 3, 3)
    end
end

@add_problem sdp function sdp_relative_entropy(handle_problem!, ::Val{test}, atol, rtol, ::Type{T}) where {T, test}
    for cplx in [false, true]
        n = 4
        X = eye(n)
        if cplx
            Y = ComplexVariable(n,n)
        else
            Y = Variable(n,n)
        end

        c1 = eye(n) in RelativeEntropyEpiCone(X, Y)
        objective = real(tr(Y))
        p = minimize(objective, c1; numeric_type = T)

        handle_problem!(p)

        if test
            @test Y.value ≈ eye(n)*exp(-1) atol=atol rtol=rtol
        end
    end
end

@add_problem sdp function sdp_geom_mean_hypocone_real_0(handle_problem!, ::Val{test}, atol, rtol, ::Type{T}) where {T, test}
    n = 4
    A = Variable(n,n)
    B = randn(n,n)
    B = B * B' # now A is positive semidefinite
    B += 0.2 * I # prevent numerical instability

    c1 = eye(n) in GeomMeanHypoCone(A, B, 0)
    objective = tr(A)
    p = minimize(objective, c1; numeric_type = T)

    handle_problem!(p)

    if test
        @test A.value ≈ eye(n) atol=atol rtol=rtol
    end
end

@add_problem sdp function sdp_geom_mean_hypocone_real_1(handle_problem!, ::Val{test}, atol, rtol, ::Type{T}) where {T, test}
    n = 4
    A = Variable(n,n)
    B = randn(n,n)
    B = B * B' # now A is positive semidefinite
    B += 0.2 * I # prevent numerical instability

    c1 = eye(n) in GeomMeanHypoCone(B, A, 1)
    objective = tr(A)
    p = minimize(objective, c1; numeric_type = T)

    handle_problem!(p)

    if test
        @test A.value ≈ eye(n) atol=atol rtol=rtol
    end
end

@add_problem sdp function sdp_geom_mean_hypocone_real_1_2(handle_problem!, ::Val{test}, atol, rtol, ::Type{T}) where {T, test}
    n = 4
    A = randn(n,n)
    A = A * A' # now A is positive semidefinite
    A += 0.2 * I # prevent numerical instability
    B = Variable(n,n)

    c1 = eye(n) in GeomMeanHypoCone(A, B, 1//2)
    objective = tr(B)
    p = minimize(objective, c1; numeric_type = T)

    handle_problem!(p)

    if test
        @test B.value ≈ A^-1 atol=atol rtol=rtol
    end
end

@add_problem sdp function sdp_geom_mean_hypocone_real_3_8(handle_problem!, ::Val{test}, atol, rtol, ::Type{T}) where {T, test}
    n = 4
    A = randn(n,n)
    A = A * A' # now A is positive semidefinite
    A += 0.2 * I # prevent numerical instability
    B = Variable(n,n)

    c1 = eye(n) in GeomMeanHypoCone(A, B, 3//8)
    objective = tr(B)
    p = minimize(objective, c1; numeric_type = T)

    handle_problem!(p)

    # A #_t B = I  =>  B = A^(1-1/t)
    if test
        @test B.value ≈ A^(-5//3) atol=atol rtol=rtol
    end
end

@add_problem sdp function sdp_geom_mean_hypocone_real_3_5(handle_problem!, ::Val{test}, atol, rtol, ::Type{T}) where {T, test}
    n = 4
    A = randn(n,n)
    A = A * A' # now A is positive semidefinite
    A += 0.2 * I # prevent numerical instability
    B = Variable(n,n)

    c1 = eye(n) in GeomMeanHypoCone(A, B, 3//5)
    objective = tr(B)
    p = minimize(objective, c1; numeric_type = T)

    handle_problem!(p)

    # A #_t B = I  =>  B = A^(1-1/t)
    if test
        @test B.value ≈ A^(-2//3) atol=atol rtol=rtol
    end
end

@add_problem sdp function sdp_geom_mean_hypocone_cplx_1_2(handle_problem!, ::Val{test}, atol, rtol, ::Type{T}) where {T, test}
    n = 4
    A = randn(ComplexF64, n, n)
    A = A * A' # now A is positive semidefinite
    A += 0.2 * I # prevent numerical instability
    B = ComplexVariable(n,n)

    c1 = eye(n) in GeomMeanHypoCone(A, B, 1//2)
    objective = real(tr(B))
    p = minimize(objective, c1; numeric_type = T)

    handle_problem!(p)

    if test
        @test real.(B.value) ≈ real.(A^-1) atol=atol rtol=rtol
        @test imag.(B.value) ≈ imag.(A^-1) atol=atol rtol=rtol
    end
end

@add_problem sdp function sdp_geom_mean_hypocone_cplx_3_8(handle_problem!, ::Val{test}, atol, rtol, ::Type{T}) where {T, test}
    n = 4
    A = randn(ComplexF64, n, n)
    A = A * A' # now A is positive semidefinite
    A += 0.2 * I # prevent numerical instability
    B = ComplexVariable(n,n)

    c1 = eye(n) in GeomMeanHypoCone(A, B, 3//8)
    objective = real(tr(B))
    p = minimize(objective, c1; numeric_type = T)

    handle_problem!(p)

    # A #_t B = I  =>  B = A^(1-1/t)
    if test
        @test real.(B.value) ≈ real.(A^(-5//3)) atol=atol rtol=rtol
        @test imag.(B.value) ≈ imag.(A^(-5//3)) atol=atol rtol=rtol
    end
end

@add_problem sdp function sdp_geom_mean_hypocone_cplx_3_5(handle_problem!, ::Val{test}, atol, rtol, ::Type{T}) where {T, test}
    n = 4
    A = randn(ComplexF64, n, n)
    A = A * A' # now A is positive semidefinite
    A += 0.2 * I # prevent numerical instability
    B = ComplexVariable(n,n)

    c1 = eye(n) in GeomMeanHypoCone(A, B, 3//5)
    objective = real(tr(B))
    p = minimize(objective, c1; numeric_type = T)

    handle_problem!(p)

    # A #_t B = I  =>  B = A^(1-1/t)
    if test
        @test real.(B.value) ≈ real.(A^(-2//3)) atol=atol rtol=rtol
        @test imag.(B.value) ≈ imag.(A^(-2//3)) atol=atol rtol=rtol
    end
end

@add_problem sdp function sdp_geom_mean_hypocone_fullhyp(handle_problem!, ::Val{test}, atol, rtol, ::Type{T}) where {T, test}
    # These matrices satisfy A^(1/2) ⪰ B and so (A,eye(2),B) \in hyp_{1/2}
    # A and B however do not satisfy A ⪰ B^2 and so [A B; B eye(2)] not
    # psd, and using fullhyp=false will give an infeasible SDP.
    A = [6.25 0; 0 16]
    B = [2 1; 1 2]
    S = Semidefinite(2)

    p = minimize(0, [B in GeomMeanHypoCone(A, eye(2), 1//2, false)]; numeric_type = T)
    handle_problem!(p)
    if test
        @test p.status == MOI.INFEASIBLE
    end

    p = minimize(0, [B in GeomMeanHypoCone(A, eye(2), 1//2)]; numeric_type = T)
    handle_problem!(p)
    if test
        @test p.status == MOI.OPTIMAL
    end
end

@add_problem sdp function sdp_geom_mean_epicone_real_neg1_optA(handle_problem!, ::Val{test}, atol, rtol, ::Type{T}) where {T, test}
    n = 3
    B = randn(n,n)
    B = B * B' # now B is positive semidefinite
    B += 0.2 * I # prevent numerical instability
    A = Variable(n,n)

    c1 = eye(n) in GeomMeanEpiCone(A, B, -1)
    objective = tr(A)
    p = maximize(objective, c1; numeric_type = T)

    handle_problem!(p)

    # A #_t B = I  =>  B = A^(1-1/t)
    if test
        @test B ≈ A.value^2 atol=atol rtol=rtol
    end
end

@add_problem sdp function sdp_geom_mean_epicone_real_neg1_optB(handle_problem!, ::Val{test}, atol, rtol, ::Type{T}) where {T, test}
    n = 3
    A = randn(n,n)
    A = A * A' # now A is positive semidefinite
    A += 0.2 * I # prevent numerical instability
    A /= tr(A) # solver has problems if B is large
    B = Variable(n,n)

    c1 = eye(n) in GeomMeanEpiCone(A, B, -1)
    objective = tr(B)
    p = minimize(objective, c1; numeric_type = T)

    handle_problem!(p)

    # A #_t B = I  =>  B = A^(1-1/t)
    if test
        @test B.value ≈ A^2 atol=atol rtol=rtol
    end
end

@add_problem sdp function sdp_geom_mean_epicone_real_neg3_5_optA(handle_problem!, ::Val{test}, atol, rtol, ::Type{T}) where {T, test}
    n = 3
    B = randn(n,n)
    B = B * B' # now B is positive semidefinite
    B += 0.2 * I # prevent numerical instability
    A = Variable(n,n)

    c1 = eye(n) in GeomMeanEpiCone(A, B, -3//5)
    objective = tr(A)
    p = maximize(objective, c1; numeric_type = T)

    handle_problem!(p)

    # A #_t B = I  =>  B = A^(1-1/t)
    if test
        @test B ≈ A.value^(8//3) atol=atol rtol=rtol
    end
end

@add_problem sdp function sdp_geom_mean_epicone_real_neg3_5_optB(handle_problem!, ::Val{test}, atol, rtol, ::Type{T}) where {T, test}
    n = 3
    A = randn(n,n)
    A = A * A' # now A is positive semidefinite
    A += 0.2 * I # prevent numerical instability
    A /= tr(A) # solver has problems if B is large
    B = Variable(n,n)

    c1 = eye(n) in GeomMeanEpiCone(A, B, -3//5)
    objective = tr(B)
    p = minimize(objective, c1; numeric_type = T)

    handle_problem!(p)

    # A #_t B = I  =>  B = A^(1-1/t)
    if test
        @test B.value ≈ A^(8//3) atol=atol rtol=rtol
    end
end

@add_problem sdp function sdp_geom_mean_epicone_real_8_5_optA(handle_problem!, ::Val{test}, atol, rtol, ::Type{T}) where {T, test}
    n = 3
    B = randn(n,n)
    B = B * B' # now B is positive semidefinite
    B += 0.2 * I # prevent numerical instability
    B /= tr(B) # solver has problems if B is large
    A = Variable(n,n)

    c1 = eye(n) in GeomMeanEpiCone(A, B, 8//5)
    objective = tr(A)
    p = minimize(objective, c1; numeric_type = T)

    handle_problem!(p)

    # A #_t B = I  =>  B = A^(1-1/t)
    if test
        @test B ≈ A.value^(3//8) atol=atol rtol=rtol
    end
end

@add_problem sdp function sdp_geom_mean_epicone_real_8_5_optB(handle_problem!, ::Val{test}, atol, rtol, ::Type{T}) where {T, test}
    n = 3
    A = randn(n,n)
    A = A * A' # now A is positive semidefinite
    A += 0.2 * I # prevent numerical instability
    B = Variable(n,n)

    c1 = eye(n) in GeomMeanEpiCone(A, B, 8//5)
    objective = tr(B)
    p = maximize(objective, c1; numeric_type = T)

    handle_problem!(p)

    # A #_t B = I  =>  B = A^(1-1/t)
    if test
        @test B.value ≈ A^(3//8) atol=atol rtol=rtol
    end
end

@add_problem sdp function sdp_quantum_relative_entropy_argcheck(handle_problem!, ::Val{test}, atol, rtol, ::Type{T}) where {T, test}
    if test
        @test_throws DimensionMismatch quantum_relative_entropy(   zeros(2,3),    zeros(2,3))
        @test_throws DimensionMismatch quantum_relative_entropy(   zeros(2,3), Variable(2,3))
        @test_throws DimensionMismatch quantum_relative_entropy(Variable(2,3),    zeros(2,3))
        @test_throws DimensionMismatch quantum_relative_entropy(Variable(2,3), Variable(2,3))

        @test_throws DimensionMismatch quantum_relative_entropy(   zeros(2,2),    zeros(3,3))
        @test_throws DimensionMismatch quantum_relative_entropy(   zeros(2,2), Variable(3,3))
        @test_throws DimensionMismatch quantum_relative_entropy(Variable(2,2),    zeros(3,3))
        @test_throws DimensionMismatch quantum_relative_entropy(Variable(2,2), Variable(2,3))

        z = zeros(2,2)
        v = Variable(2,2)
        nh = [ 1 1 ; 0 1 ] # not hermitian
        np = [ 1 2 ; 2 1 ] # not positive semidefinite
        @test_throws DomainError quantum_relative_entropy(z, nh)
        @test_throws DomainError quantum_relative_entropy(z, np)
        @test_throws DomainError quantum_relative_entropy(nh, z)
        @test_throws DomainError quantum_relative_entropy(np, z)
        @test_throws DomainError quantum_relative_entropy(v, nh)
        @test_throws DomainError quantum_relative_entropy(v, np)
        @test_throws DomainError quantum_relative_entropy(nh, v)
        @test_throws DomainError quantum_relative_entropy(np, v)
    end
end

function sdp_quantum_relative_entropy_impl(handle_problem!, ::Val{test}, atol, rtol, ::Type{T}, lowrank::Bool, mode::Integer) where {T, test}
    # Too slow for n>3.  Anything that can be done to speed it up?
    n = 3
    X = randn(ComplexF64, n, n)
    if lowrank
        X[:,1] .= 0
    end
    X = X * X' # now X is positive semidefinite
    X /= tr(X) # make it a quantum state

    A = ComplexVariable(n,n)
    B = ComplexVariable(n,n)
    constraints = Array{Constraint, 1}([
        A == X,
    ])
    if mode != 5
        push!(constraints, tr(B) == 1)
    end

    if mode == 1
        objective = quantum_relative_entropy(A, B)
    elseif mode == 2
        objective = quantum_relative_entropy(X, B)
    elseif mode == 3
        objective = quantum_relative_entropy(B, A)
    elseif mode == 4
        objective = quantum_relative_entropy(B, X)
    elseif mode == 5
        objective = quantum_relative_entropy(X, X)
    end

    p = minimize(objective, constraints; numeric_type = T)
    handle_problem!(p)

    if test
        @test real.(evaluate(A)) ≈ real.(X) atol=atol rtol=rtol
        @test imag.(evaluate(A)) ≈ imag.(X) atol=atol rtol=rtol
        if mode != 5
            @test real.(evaluate(B)) ≈ real.(X) atol=atol rtol=rtol
            @test imag.(evaluate(B)) ≈ imag.(X) atol=atol rtol=rtol
        end
        @test p.optval ≈ 0 atol=atol rtol=rtol
        if mode == 1
            @test p.optval ≈ evaluate(quantum_relative_entropy(A, B)) atol=atol rtol=rtol
        elseif mode == 2
            @test p.optval ≈ evaluate(quantum_relative_entropy(X, B)) atol=atol rtol=rtol
        elseif mode == 3
            @test p.optval ≈ evaluate(quantum_relative_entropy(B, A)) atol=atol rtol=rtol
        elseif mode == 4
            @test p.optval ≈ evaluate(quantum_relative_entropy(B, X)) atol=atol rtol=rtol
        elseif mode == 5
            @test p.optval ≈ evaluate(quantum_relative_entropy(X, X)) atol=atol rtol=rtol
        end
    end
end

@add_problem sdp function sdp_quantum_relative_entropy1_fullrank(handle_problem!, ::Val{test}, atol, rtol, ::Type{T}) where {T, test}
    sdp_quantum_relative_entropy_impl(handle_problem!, Val(test), atol, rtol, T, false, 1)
end

@add_problem sdp function sdp_quantum_relative_entropy2_fullrank(handle_problem!, ::Val{test}, atol, rtol, ::Type{T}) where {T, test}
    sdp_quantum_relative_entropy_impl(handle_problem!, Val(test), atol, rtol, T, false, 2)
end

@add_problem sdp function sdp_quantum_relative_entropy3_fullrank(handle_problem!, ::Val{test}, atol, rtol, ::Type{T}) where {T, test}
    sdp_quantum_relative_entropy_impl(handle_problem!, Val(test), atol, rtol, T, false, 3)
end

@add_problem sdp function sdp_quantum_relative_entropy4_fullrank(handle_problem!, ::Val{test}, atol, rtol, ::Type{T}) where {T, test}
    sdp_quantum_relative_entropy_impl(handle_problem!, Val(test), atol, rtol, T, false, 4)
end

@add_problem sdp function sdp_quantum_relative_entropy5_fullrank(handle_problem!, ::Val{test}, atol, rtol, ::Type{T}) where {T, test}
    sdp_quantum_relative_entropy_impl(handle_problem!, Val(test), atol, rtol, T, false, 5)
end

@add_problem sdp function sdp_quantum_relative_entropy1_lowrank(handle_problem!, ::Val{test}, atol, rtol, ::Type{T}) where {T, test}
    sdp_quantum_relative_entropy_impl(handle_problem!, Val(test), atol, rtol, T, true, 1)
end

@add_problem sdp function sdp_quantum_relative_entropy2_lowrank(handle_problem!, ::Val{test}, atol, rtol, ::Type{T}) where {T, test}
    sdp_quantum_relative_entropy_impl(handle_problem!, Val(test), atol, rtol, T, true, 2)
end

@add_problem sdp function sdp_quantum_relative_entropy3_lowrank(handle_problem!, ::Val{test}, atol, rtol, ::Type{T}) where {T, test}
    sdp_quantum_relative_entropy_impl(handle_problem!, Val(test), atol, rtol, T, true, 3)
end

@add_problem sdp function sdp_quantum_relative_entropy4_lowrank(handle_problem!, ::Val{test}, atol, rtol, ::Type{T}) where {T, test}
    sdp_quantum_relative_entropy_impl(handle_problem!, Val(test), atol, rtol, T, true, 4)
end

@add_problem sdp function sdp_quantum_relative_entropy5_lowrank(handle_problem!, ::Val{test}, atol, rtol, ::Type{T}) where {T, test}
    sdp_quantum_relative_entropy_impl(handle_problem!, Val(test), atol, rtol, T, true, 5)
end

@add_problem sdp function sdp_quantum_relative_entropy_const(handle_problem!, ::Val{test}, atol, rtol, ::Type{T}) where {T, test}
    if test
        @test quantum_relative_entropy(diagm(0=>[0.5, 0.5]), diagm(0=>[0.5, 0.5])) ≈ 0 atol=atol rtol=rtol
        @test quantum_relative_entropy(diagm(0=>[1.0, 0.0]), diagm(0=>[1.0, 0.0])) ≈ 0 atol=atol rtol=rtol
        @test quantum_relative_entropy(diagm(0=>[1.0, 0.0]), diagm(0=>[0.5, 0.5])) ≈ log(2) atol=atol rtol=rtol
        @test isinf(quantum_relative_entropy(diagm(0=>[0.5, 0.5]), diagm(0=>[1.0, 0.0])))
    end
end

@add_problem sdp function sdp_quantum_channel_capacity(handle_problem!, ::Val{test}, atol, rtol, ::Type{T}) where {T, test}
    # Capacity of classical-quantum channel.  See chapter 3 of arxiv:1705.06671,
    # Fawzi, Fawzi "Efficient optimization of the quantum relative entropy".
    # For channels with two symbols, both pure states, we have an exact formula to compare against.

    n = 3

    for cplx in [false, true]
        if cplx
            ψ1 = normalize(randn(ComplexF64, n))
            ψ2 = normalize(randn(ComplexF64, n))
        else
            ψ1 = normalize(randn(Float64, n))
            ψ2 = normalize(randn(Float64, n))
        end
        ρ1 = ψ1 * ψ1'
        ρ2 = ψ2 * ψ2'

        p1 = Variable()
        p2 = Variable()

        objective = (quantum_entropy(p1*ρ1 + p2*ρ2) - p1*quantum_entropy(ρ1) - p2*quantum_entropy(ρ2)) / log(2)
        constraints = [ p1 >= 0, p2 >= 0, p1+p2 == 1 ]
        p = maximize(objective, constraints; numeric_type = T)
        handle_problem!(p)

        ϵ = abs(ψ1' * ψ2)
        q = (1 + ϵ)/2
        r = [q, 1-q]
        v = -r' * log2.(r)

        if test
            @test p.optval ≈ v atol=atol rtol=rtol
            @test evaluate(objective) ≈ v atol=atol rtol=rtol

            @test_throws DimensionMismatch quantum_entropy(Variable(2, 3))
        end
    end
end

@add_problem sdp function sdp_trace_mpower_argcheck(handle_problem!, ::Val{test}, atol, rtol, ::Type{T}) where {T, test}
    if test
        @test_throws DimensionMismatch trace_mpower(   zeros(2,3), 1//2, zeros(2,3))
        @test_throws DimensionMismatch trace_mpower(Variable(2,3), 1//2, zeros(2,3))

        @test_throws DimensionMismatch trace_mpower(   zeros(2,2), 1//2, zeros(3,3))
        @test_throws DimensionMismatch trace_mpower(Variable(2,2), 1//2, zeros(3,3))

        @test_throws DomainError trace_mpower(Variable(3,3),  5//2, zeros(3,3))

        nh = [ 1 1 ; 0 1 ] # not hermitian
        np = [ 1 2 ; 2 1 ] # not positive semidefinite
        @test_throws DomainError trace_mpower(Variable(2,2), 1//2, nh)
        @test_throws DomainError trace_mpower(Variable(2,2), 1//2, np)
    end
end

@add_problem sdp function sdp_trace_mpower_real_2_3(handle_problem!, ::Val{test}, atol, rtol, ::Type{T}) where {T, test}
    # Too slow for n>3.  Anything that can be done to speed it up?
    n = 3
    A = randn(n,n)
    A = A * A' # now A is positive semidefinite
    C = randn(n,n)
    C = C * C' # now C is positive semidefinite
    t = 2//3

    B = Semidefinite(n,n)

    constraints = [
        B ⪯ A
    ]
    objective = trace_mpower(B, t, C)
    p = maximize(objective, constraints; numeric_type = T)

    handle_problem!(p)

    if test
        @test real.(B.value) ≈ real.(A) atol=atol rtol=rtol
        @test imag.(B.value) ≈ imag.(A) atol=atol rtol=rtol
        @test p.optval ≈ tr(C*A^t) atol=atol rtol=rtol
        @test p.optval ≈ trace_mpower(A, t, C) atol=atol rtol=rtol
        @test p.optval ≈ evaluate(trace_mpower(B, t, C)) atol=atol rtol=rtol
    end
end

@add_problem sdp function sdp_trace_mpower_cplx_2_3(handle_problem!, ::Val{test}, atol, rtol, ::Type{T}) where {T, test}
    # Too slow for n>3.  Anything that can be done to speed it up?
    n = 3
    A = randn(ComplexF64, n, n)
    A = A * A' # now A is positive semidefinite
    C = randn(ComplexF64, n, n)
    C = C * C' # now C is positive semidefinite
    t = 2//3

    B = HermitianSemidefinite(n,n)

    constraints = [
        B ⪯ A
    ]
    objective = trace_mpower(B, t, C)
    p = maximize(objective, constraints; numeric_type = T)

    handle_problem!(p)

    if test
        @test real.(B.value) ≈ real.(A) atol=atol rtol=rtol
        @test imag.(B.value) ≈ imag.(A) atol=atol rtol=rtol
        @test p.optval ≈ tr(C*A^t) atol=atol rtol=rtol
        @test p.optval ≈ trace_mpower(A, t, C) atol=atol rtol=rtol
        @test p.optval ≈ evaluate(trace_mpower(B, t, C)) atol=atol rtol=rtol
    end
end

@add_problem sdp function sdp_trace_mpower_real_5_4(handle_problem!, ::Val{test}, atol, rtol, ::Type{T}) where {T, test}
    # Too slow for n>3.  Anything that can be done to speed it up?
    n = 3
    A = randn(n,n)
    A = A * A' # now A is positive semidefinite
    C = randn(n,n)
    C = C * C' # now C is positive semidefinite
    t = 5//4

    B = Semidefinite(n,n)

    constraints = [
        # note: B ⪰ A doesn't work for this unit test because trace_mpower is not operator monotone
        # unless t ∈ [0,1].
        B == A
    ]
    objective = trace_mpower(B, t, C)
    p = minimize(objective, constraints; numeric_type = T)

    handle_problem!(p)

    if test
        @test real.(B.value) ≈ real.(A) atol=atol rtol=rtol
        @test imag.(B.value) ≈ imag.(A) atol=atol rtol=rtol
        @test p.optval ≈ tr(C*A^t) atol=atol rtol=rtol
        @test p.optval ≈ trace_mpower(A, t, C) atol=atol rtol=rtol
        @test p.optval ≈ evaluate(trace_mpower(B, t, C)) atol=atol rtol=rtol
    end
end

@add_problem sdp function sdp_trace_mpower_cplx_5_4(handle_problem!, ::Val{test}, atol, rtol, ::Type{T}) where {T, test}
    # Too slow for n>3.  Anything that can be done to speed it up?
    n = 3
    A = randn(ComplexF64, n, n)
    A = A * A' # now A is positive semidefinite
    C = randn(ComplexF64, n, n)
    C = C * C' # now C is positive semidefinite
    t = 5//4

    B = HermitianSemidefinite(n,n)

    constraints = [
        # note: B ⪰ A doesn't work for this unit test because trace_mpower is not operator monotone
        # unless t ∈ [0,1].
        B == A
    ]
    objective = trace_mpower(B, t, C)
    p = minimize(objective, constraints; numeric_type = T)

    handle_problem!(p)

    if test
        @test real.(B.value) ≈ real.(A) atol=atol rtol=rtol
        @test imag.(B.value) ≈ imag.(A) atol=atol rtol=rtol
        @test p.optval ≈ tr(C*A^t) atol=atol rtol=rtol
        @test p.optval ≈ trace_mpower(A, t, C) atol=atol rtol=rtol
        @test p.optval ≈ evaluate(trace_mpower(B, t, C)) atol=atol rtol=rtol
    end
end

@add_problem sdp function sdp_trace_mpower_real_neg1_4(handle_problem!, ::Val{test}, atol, rtol, ::Type{T}) where {T, test}
    # Too slow for n>3.  Anything that can be done to speed it up?
    n = 3
    A = randn(n,n)
    A = A * A' # now A is positive semidefinite
    A += 0.2 * I # prevent numerical instability
    C = randn(n,n)
    C = C * C' # now C is positive semidefinite
    t = -1//4

    B = Semidefinite(n,n)

    constraints = [
        # -trace_mpower is operator monotone for t ∈ [-1, 0].
        B ⪯ A
    ]
    objective = trace_mpower(B, t, C)
    p = minimize(objective, constraints; numeric_type = T)

    handle_problem!(p)

    if test
        @test real.(B.value) ≈ real.(A) atol=atol rtol=rtol
        @test imag.(B.value) ≈ imag.(A) atol=atol rtol=rtol
        @test p.optval ≈ tr(C*A^t) atol=atol rtol=rtol
        @test p.optval ≈ trace_mpower(A, t, C) atol=atol rtol=rtol
        @test p.optval ≈ evaluate(trace_mpower(B, t, C)) atol=atol rtol=rtol
    end
end

@add_problem sdp function sdp_trace_mpower_cplx_neg1_4(handle_problem!, ::Val{test}, atol, rtol, ::Type{T}) where {T, test}
    # Too slow for n>3.  Anything that can be done to speed it up?
    n = 3
    A = randn(ComplexF64, n, n)
    A = A * A' # now A is positive semidefinite
    A += 0.2 * I # prevent numerical instability
    C = randn(ComplexF64, n, n)
    C = C * C' # now C is positive semidefinite
    t = -1//4

    B = HermitianSemidefinite(n,n)

    constraints = [
        # -trace_mpower is operator monotone for t ∈ [-1, 0].
        B ⪯ A
    ]
    objective = trace_mpower(B, t, C)
    p = minimize(objective, constraints; numeric_type = T)

    handle_problem!(p)

    if test
        @test real.(B.value) ≈ real.(A) atol=atol rtol=rtol
        @test imag.(B.value) ≈ imag.(A) atol=atol rtol=rtol
        @test p.optval ≈ tr(C*A^t) atol=atol rtol=rtol
        @test p.optval ≈ trace_mpower(A, t, C) atol=atol rtol=rtol
        @test p.optval ≈ evaluate(trace_mpower(B, t, C)) atol=atol rtol=rtol
    end
end

@add_problem sdp function sdp_trace_logm_argcheck(handle_problem!, ::Val{test}, atol, rtol, ::Type{T}) where {T, test}
    if test
        @test_throws DimensionMismatch trace_logm(   zeros(2,3), zeros(2,3))
        @test_throws DimensionMismatch trace_logm(Variable(2,3), zeros(2,3))

        @test_throws DimensionMismatch trace_logm(   zeros(2,2), zeros(3,3))
        @test_throws DimensionMismatch trace_logm(Variable(2,2), zeros(3,3))

        nh = [ 1 1 ; 0 1 ] # not hermitian
        np = [ 1 2 ; 2 1 ] # not positive semidefinite
        @test_throws DomainError trace_logm(Variable(2,2), nh)
        @test_throws DomainError trace_logm(Variable(2,2), np)
    end
end

@add_problem sdp function sdp_trace_logm_real(handle_problem!, ::Val{test}, atol, rtol, ::Type{T}) where {T, test}
    n = 3
    C = randn(Float64, n, n)
    C = C * C' # now C is positive semidefinite

    X = Variable(n,n)

    constraints = [
        X ⪯ eye(n)
    ]
    objective = trace_logm(X, C)
    p = maximize(objective, constraints; numeric_type = T)

    handle_problem!(p)

    if test
        @test real.(X.value) ≈ real.(eye(n)) atol=atol rtol=rtol
        @test imag.(X.value) ≈ imag.(eye(n)) atol=atol rtol=rtol
        @test p.optval ≈ tr(C*log(evaluate(X))) atol=atol rtol=rtol
        @test p.optval ≈ trace_logm(evaluate(X), C) atol=atol rtol=rtol
        @test p.optval ≈ evaluate(trace_logm(X, C)) atol=atol rtol=rtol
    end
end

@add_problem sdp function sdp_trace_logm_cplx(handle_problem!, ::Val{test}, atol, rtol, ::Type{T}) where {T, test}
    n = 3
    C = randn(ComplexF64, n, n)
    C = C * C' # now C is positive semidefinite

    X = ComplexVariable(n,n)

    constraints = [
        X ⪯ eye(n)
    ]
    objective = trace_logm(X, C)
    p = maximize(objective, constraints; numeric_type = T)

    handle_problem!(p)

    if test
        @test real.(X.value) ≈ real.(eye(n)) atol=atol rtol=rtol
        @test imag.(X.value) ≈ imag.(eye(n)) atol=atol rtol=rtol
        @test p.optval ≈ tr(C*log(evaluate(X))) atol=atol rtol=rtol
        @test p.optval ≈ trace_logm(evaluate(X), C) atol=atol rtol=rtol
        @test p.optval ≈ evaluate(trace_logm(X, C)) atol=atol rtol=rtol
    end
end

@add_problem sdp function sdp_lieb_ando(handle_problem!, ::Val{test}, atol, rtol, ::Type{T}) where {T, test}
    if test
        @test_throws DomainError lieb_ando(   zeros(3,3),    zeros(3,3), eye(3), -3//2)
        @test_throws DomainError lieb_ando(   zeros(3,3),    zeros(3,3), eye(3), 5//2)
        @test_throws DomainError lieb_ando(   zeros(3,3), Variable(3,3), eye(3), -3//2)
        @test_throws DomainError lieb_ando(   zeros(3,3), Variable(3,3), eye(3), 5//2)
        @test_throws DomainError lieb_ando(Variable(3,3),    zeros(3,3), eye(3), -3//2)
        @test_throws DomainError lieb_ando(Variable(3,3),    zeros(3,3), eye(3), 5//2)
        @test_throws DomainError lieb_ando(Variable(3,3), Variable(3,3), eye(3), -3//2)
        @test_throws DomainError lieb_ando(Variable(3,3), Variable(3,3), eye(3), 5//2)
    end

    for n in [2, 3]
        for t in [1//2, 1//4, 3//4, 1//8, 3//2, 5//4]
            for cplx in [false, true]
                #@show n,t,cplx

                if cplx
                    A = randn(ComplexF64, n, n)
                    B = randn(ComplexF64, n, n)
                    X = ComplexVariable(n, n)
                    Y = ComplexVariable(n, n)
                else
                    A = randn(n, n)
                    B = randn(n, n)
                    X = Variable(n, n)
                    Y = Variable(n, n)
                end
                A = A*A'
                B = B*B'
                A = A/tr(A)
                B = B/tr(B)

                QtAB = lieb_ando(A,B,eye(n),t)

                objective = lieb_ando(X,B,eye(n),t)
                if t >= 0 && t <= 1
                    p = maximize(objective, [X == A]; numeric_type = T)
                else
                    p = minimize(objective, [X == A]; numeric_type = T)
                end
                handle_problem!(p)
                if test
                    @test p.optval ≈ QtAB atol=atol*5 rtol=rtol
                end

                objective = lieb_ando(A,Y,eye(n),t)
                if t >= 0 && t <= 1
                    p = maximize(objective, [Y == B]; numeric_type = T)
                else
                    p = minimize(objective, [Y == B]; numeric_type = T)
                end
                handle_problem!(p)
                if test
                    @test p.optval ≈ QtAB atol=atol*5 rtol=rtol
                end

                objective = lieb_ando(X,Y,eye(n),t)
                if t >= 0 && t <= 1
                    p = maximize(objective, [X == A, Y == B]; numeric_type = T)
                else
                    p = minimize(objective, [X == A, Y == B]; numeric_type = T)
                end
                handle_problem!(p)
                if test
                    @test p.optval ≈ QtAB atol=atol*5 rtol=rtol
                end
            end
        end
    end
end

@add_problem sdp function sdp_min_maxeig_canon_lmi(handle_problem!, ::Val{test}, atol, rtol, ::Type{T}) where {T, test}
    # Minimize the maximum eigenvalue of an affine matrix function of a vector argument. Formulated
    # as a linear optimization subject to a canonical LMI (SDP) constraint:
    #
    # minimize      λ
    # subject to    A(x)⪯λ
    # where x is a vector in Rⁿ and the matrix function A(x) = A₀ + A₁x₁ + A₂x₂ + … + Aₙxₙ.
    #
    # Besides serving as a benchmark, it also tests the solution of the Issue
    # https://github.com/jump-dev/Convex.jl/issues/447
    # (the issues boils down to  evaluation of `A₁x₁`, because `x₁` was a matrix of size (1,1)).

    A₀ = [1.0 2.0; 2.0 3.0]
    A₁ = [1.0 0.0; 0.0 -1.0]
    A₂ = [4.0 5.0; 5.0 -6.0]

    x = Variable(2)
    λ = Variable()

    A = A₀ + A₁*x[1] + A₂*x[2]

    p = minimize(λ,A⪯λ)

    handle_problem!(p)

    if test
        @test λ.value ≈ 2.0 atol=atol rtol=rtol
        @test x.value ≈ [1.0, 0.0] atol=atol rtol=rtol
        @test evaluate(A) ≈ [2.0 2.0; 2.0 2.0] atol=atol rtol=rtol
    end
end
