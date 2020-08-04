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
    p = minimize(sumlargesteigs(x, 2), [x[i,:] >= i for i=1:3]...; numeric_type = T)

    handle_problem!(p)

    if test
        @test p.optval ≈ 8.4853 atol=atol rtol=rtol
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
        @test real(x1) ≈ evaluate(xr) atol=atol rtol=rtol
        @test imag(x1) ≈ evaluate(xi) atol=atol rtol=rtol
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

        @test real(evaluate(x)) ≈ real(a) atol=atol rtol=rtol
        @test imag(evaluate(x)) ≈ imag(a) atol=atol rtol=rtol
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

        @test real.(evaluate(x)) ≈ real.(a) atol=atol rtol=rtol
        @test imag.(evaluate(x)) ≈ imag.(a) atol=atol rtol=rtol
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

        @test real(evaluate(x)) ≈ real(a) atol=atol rtol=rtol
        @test imag(evaluate(x)) ≈ imag(a) atol=atol rtol=rtol
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
        @test real.(evaluate(x)) ≈ real.(posA) atol=atol rtol=rtol
        @test imag.(evaluate(x)) ≈ imag.(posA) atol=atol rtol=rtol
    end
end
