# Copyright (c) 2014: Madeleine Udell and contributors
#
# Use of this source code is governed by a BSD-style license that can be found
# in the LICENSE file or at https://opensource.org/license/bsd-2-clause

# TODO: uncomment vexity checks once SDP on vars/constraints changes vexity of problem
@add_problem sdp function sdp_sdp_variables(
    handle_problem!,
    ::Val{test},
    atol,
    rtol,
    ::Type{T},
) where {T,test}
    y = Semidefinite((2, 2))
    p = minimize(y[1, 1]; numeric_type = T)

    # @fact problem_vexity(p) --> ConvexVexity()
    handle_problem!(p)
    if test
        @test p.optval ≈ 0 atol = atol rtol = rtol
    end

    y = Semidefinite((3, 3))
    p = minimize(y[1, 1], y[2, 2] == 1; numeric_type = T)

    # @fact problem_vexity(p) --> ConvexVexity()
    handle_problem!(p)
    if test
        @test p.optval ≈ 0 atol = atol rtol = rtol
    end

    # Solution is obtained as y[2,2] -> infinity
    # This test fails on Mosek. See
    # https://github.com/JuliaOpt/Mosek.jl/issues/29
    # y = Semidefinite((2, 2))
    # p = minimize(y[1, 1], y[1, 2] == 1; numeric_type = T)

    # # @fact problem_vexity(p) --> ConvexVexity()
    # handle_problem!(p)
    # @fact p.optval --> roughly(0, atol)

    y = Semidefinite(3)
    p = minimize(sum(diag(y)), y[1, 1] == 1; numeric_type = T)

    # @fact problem_vexity(p) --> ConvexVexity()
    handle_problem!(p)
    if test
        @test p.optval ≈ 1 atol = atol rtol = rtol
    end

    y = Semidefinite((3, 3))
    p = minimize(tr(y), y[2, 1] <= 4, y[2, 2] >= 3; numeric_type = T)

    # @fact problem_vexity(p) --> ConvexVexity()
    handle_problem!(p)
    if test
        @test p.optval ≈ 3 atol = atol rtol = rtol
    end

    x = Variable(Positive())
    y = Semidefinite(3)
    p = minimize(y[1, 2], y[2, 1] == 1; numeric_type = T)

    # @fact problem_vexity(p) --> ConvexVexity()
    handle_problem!(p)
    if test
        @test p.optval ≈ 1 atol = atol rtol = rtol
    end
end

@add_problem sdp function sdp_sdp_constraints(
    handle_problem!,
    ::Val{test},
    atol,
    rtol,
    ::Type{T},
) where {T,test}
    x = Variable(Positive())
    y = Variable((3, 3))
    p = minimize(x + y[1, 1], y ⪰ 0, x >= 1, y[2, 1] == 1; numeric_type = T)
    # @fact problem_vexity(p) --> ConvexVexity()
    handle_problem!(p)
    if test
        @test p.optval ≈ 1 atol = atol rtol = rtol
    end
end

@add_problem sdp function sdp_nuclear_norm_atom(
    handle_problem!,
    ::Val{test},
    atol,
    rtol,
    ::Type{T},
) where {T,test}
    y = Semidefinite(3)
    p = minimize(
        nuclearnorm(y),
        y[2, 1] <= 4,
        y[2, 2] >= 3,
        y[3, 3] <= 2;
        numeric_type = T,
    )

    if test
        @test problem_vexity(p) == ConvexVexity()
    end
    handle_problem!(p)
    if test
        @test p.optval ≈ 3 atol = atol rtol = rtol
        @test evaluate(nuclearnorm(y)) ≈ 3 atol = atol rtol = rtol
    end
end

@add_problem sdp function sdp_nuclear_norm_atom_complex(
    handle_problem!,
    ::Val{test},
    atol,
    rtol,
    ::Type{T},
) where {T,test}
    A = [1 2im 3 4; 4im 3im 2 1; 4 5 6 7]
    y = ComplexVariable(3, 4)
    p = minimize(nuclearnorm(y), y == A; numeric_type = T)

    if test
        @test problem_vexity(p) == ConvexVexity()
    end
    handle_problem!(p)
    if test
        @test p.optval ≈ sum(LinearAlgebra.svdvals(A)) atol = atol rtol = rtol
        @test evaluate(nuclearnorm(y)) ≈ sum(LinearAlgebra.svdvals(A)) atol =
            atol rtol = rtol
    end
end

@add_problem sdp function sdp_operator_norm_atom(
    handle_problem!,
    ::Val{test},
    atol,
    rtol,
    ::Type{T},
) where {T,test}
    y = Variable((3, 3))
    p = minimize(
        LinearAlgebra.opnorm(y),
        y[2, 1] <= 4,
        y[2, 2] >= 3,
        sum(y) >= 12;
        numeric_type = T,
    )

    if test
        @test problem_vexity(p) == ConvexVexity()
    end
    handle_problem!(p)
    if test
        @test p.optval ≈ 4 atol = atol rtol = rtol
        @test evaluate(LinearAlgebra.opnorm(y)) ≈ 4 atol = atol rtol = rtol
    end
end

@add_problem sdp function sdp_operator_norm_atom_complex(
    handle_problem!,
    ::Val{test},
    atol,
    rtol,
    ::Type{T},
) where {T,test}
    A = [1 2im 3 4; 4im 3im 2 1; 4 5 6 7]
    y = ComplexVariable(3, 4)
    p = minimize(LinearAlgebra.opnorm(y), y == A; numeric_type = T)

    if test
        @test problem_vexity(p) == ConvexVexity()
    end
    handle_problem!(p)
    if test
        @test p.optval ≈ maximum(LinearAlgebra.svdvals(A)) atol = atol rtol =
            rtol
        @test evaluate(LinearAlgebra.opnorm(y)) ≈
              maximum(LinearAlgebra.svdvals(A)) atol = atol rtol = rtol
    end
end

@add_problem sdp function sdp_sigma_max_atom(
    handle_problem!,
    ::Val{test},
    atol,
    rtol,
    ::Type{T},
) where {T,test}
    y = Variable((3, 3))
    p = minimize(
        sigmamax(y),
        y[2, 1] <= 4,
        y[2, 2] >= 3,
        sum(y) >= 12;
        numeric_type = T,
    )

    if test
        @test problem_vexity(p) == ConvexVexity()
    end
    handle_problem!(p)
    if test
        @test p.optval ≈ 4 atol = atol rtol = rtol
        @test evaluate(sigmamax(y)) ≈ 4 atol = atol rtol = rtol
    end
end

@add_problem sdp function sdp_dual_lambda_max_atom(
    handle_problem!,
    ::Val{test},
    atol,
    rtol,
    ::Type{T},
) where {T,test}
    y = Semidefinite(3)
    p = minimize(eigmax(y), y[1, 1] >= 4; numeric_type = T)

    if test
        @test problem_vexity(p) == ConvexVexity()
    end
    handle_problem!(p)
    if test
        @test p.optval ≈ 4 atol = atol rtol = rtol
        @test evaluate(eigmax(y)) ≈ 4 atol = atol rtol = rtol
    end

    # https://github.com/jump-dev/Convex.jl/issues/337
    x = ComplexVariable(2, 2)
    p = minimize(
        eigmax(x),
        [x[1, 2] == im, x[2, 2] == 1.0, x ⪰ -eye(2)];
        numeric_type = T,
    )
    handle_problem!(p)
    if test
        @test p.optval ≈ 1.5 atol = atol rtol = rtol
        @test p.constraints[1].dual ≈ im atol = atol rtol = rtol
        @test p.constraints[2].dual ≈ 0.75 atol = atol rtol = rtol
    end
end

@add_problem sdp function sdp_lambda_min_atom(
    handle_problem!,
    ::Val{test},
    atol,
    rtol,
    ::Type{T},
) where {T,test}
    y = Semidefinite(3)
    p = maximize(eigmin(y), tr(y) <= 6; numeric_type = T)

    if test
        @test problem_vexity(p) == ConvexVexity()
    end
    handle_problem!(p)
    if test
        @test p.optval ≈ 2 atol = atol rtol = rtol
        @test evaluate(eigmin(y)) ≈ 2 atol = atol rtol = rtol
    end
end

@add_problem sdp function sdp_matrix_frac_atom(
    handle_problem!,
    ::Val{test},
    atol,
    rtol,
    ::Type{T},
) where {T,test}
    x = [1, 2, 3]
    P = Variable(3, 3)
    p = minimize(
        matrixfrac(x, P),
        P <= 2 * eye(3),
        P >= 0.5 * eye(3);
        numeric_type = T,
    )

    if test
        @test problem_vexity(p) == ConvexVexity()
    end
    handle_problem!(p)
    if test
        @test p.optval ≈ 7 atol = atol rtol = rtol
        @test (evaluate(matrixfrac(x, P)))[1] ≈ 7 atol = atol rtol = rtol
    end
end

@add_problem sdp function sdp_matrix_frac_atom_both_arguments_variable(
    handle_problem!,
    ::Val{test},
    atol,
    rtol,
    ::Type{T},
) where {T,test}
    x = Variable(3)
    P = Variable(3, 3)
    p = minimize(matrixfrac(x, P), eigmax(P) <= 2, x[1] >= 1; numeric_type = T)

    if test
        @test problem_vexity(p) == ConvexVexity()
    end
    handle_problem!(p)
    if test
        @test p.optval ≈ 0.5 atol = atol rtol = rtol
        @test (evaluate(matrixfrac(x, P)))[1] ≈ 0.5 atol = atol rtol = rtol
    end
end

@add_problem sdp function sdp_sum_largest_eigs(
    handle_problem!,
    ::Val{test},
    atol,
    rtol,
    ::Type{T},
) where {T,test}
    x = Semidefinite(3)
    p = minimize(sumlargesteigs(x, 2), x >= 1; numeric_type = T)

    handle_problem!(p)

    if test
        @test p.optval ≈ 3 atol = atol rtol = rtol
        @test evaluate(x) ≈ ones(3, 3) atol = atol rtol = rtol
    end

    x = Semidefinite(3)
    p = minimize(
        sumlargesteigs(x, 2) + sumlargesteigs(x, 0),
        [x[i, :] >= i for i in 1:3]...;
        numeric_type = T,
    )

    handle_problem!(p)

    if test
        @test p.optval ≈ 8.4853 atol = atol rtol = rtol
    end

    A = [
        1 -2im 3 4
        2im 7 3im 5
        3 -3im 2 9
        4 5 9 4
    ]

    x = ComplexVariable(4, 4)
    p = minimize(sumlargesteigs(x, 3), x == A; numeric_type = T)

    handle_problem!(p)

    if test
        @test p.optval ≈ sum(LinearAlgebra.eigvals(A)[2:end]) atol = atol rtol =
            rtol
    end

    x1 = Semidefinite(3)
    p1 = minimize(eigmax(x1), x1[1, 1] >= 4; numeric_type = T)

    handle_problem!(p1)

    x2 = Semidefinite(3)
    p2 = minimize(sumlargesteigs(x2, 1), x2[1, 1] >= 4; numeric_type = T)

    handle_problem!(p2)

    if test
        @test p1.optval ≈ p2.optval atol = atol rtol = rtol
    end

    x1 = Semidefinite(3)
    p1 = minimize(eigmax(x1), [x1[i, :] >= i for i in 1:3]...; numeric_type = T)

    handle_problem!(p1)

    x2 = Semidefinite(3)
    p2 = minimize(
        sumlargesteigs(x2, 1),
        [x2[i, :] >= i for i in 1:3]...;
        numeric_type = T,
    )

    handle_problem!(p2)

    if test
        @test p1.optval ≈ p2.optval atol = atol rtol = rtol
    end
end

@add_problem sdp function sdp_kron_atom(
    handle_problem!,
    ::Val{test},
    atol,
    rtol,
    ::Type{T},
) where {T,test}
    id = eye(4)
    X = Semidefinite(4)
    W = kron(id, X)
    p = maximize(tr(W), tr(X) ≤ 1; numeric_type = T)

    if test
        @test problem_vexity(p) == AffineVexity()
    end
    handle_problem!(p)
    if test
        @test p.optval ≈ 4 atol = atol rtol = rtol
    end
end

@add_problem sdp function sdp_Partial_trace(
    handle_problem!,
    ::Val{test},
    atol,
    rtol,
    ::Type{T},
) where {T,test}
    A = Semidefinite(2)
    B = [1 0; 0 0]
    ρ = kron(B, A)
    constraints = [
        partialtrace(ρ, 1, [2; 2]) ==
        [0.09942819 0.29923607; 0.29923607 0.90057181],
        isposdef(ρ),
    ]
    p = satisfy(constraints; numeric_type = T)

    handle_problem!(p)
    if test
        @test evaluate(ρ) ≈ [
            0.09942819 0.29923607 0 0
            0.299237 0.900572 0 0
            0 0 0 0
            0 0 0 0
        ] atol = atol rtol = rtol
        @test evaluate(partialtrace(ρ, 1, [2; 2])) ≈
              [0.09942819 0.29923607; 0.29923607 0.90057181] atol = atol rtol =
            rtol
    end

    function rand_normalized(n)
        A = 5 * randn(n, n) + im * 5 * randn(n, n)
        return A / tr(A)
    end

    As = [rand_normalized(3) for _ in 1:5]
    Bs = [rand_normalized(2) for _ in 1:5]
    p = rand(5)

    AB = sum(i -> p[i] * kron(As[i], Bs[i]), 1:5)
    if test
        @test partialtrace(AB, 2, [3, 2]) ≈ sum(p .* As) atol = atol rtol = rtol
        @test partialtrace(AB, 1, [3, 2]) ≈ sum(p .* Bs) atol = atol rtol = rtol
    end

    A, B, C = rand(5, 5), rand(4, 4), rand(3, 3)
    ABC = kron(kron(A, B), C)
    if test
        @test kron(A, B) * tr(C) ≈ partialtrace(ABC, 3, [5, 4, 3]) atol = atol rtol =
            rtol
    end

    # Test 281
    A = rand(6, 6)
    expr = partialtrace(constant(A), 1, [2, 3])
    if test
        @test size(expr) == size(evaluate(expr))

        @test_throws ArgumentError partialtrace(rand(6, 6), 3, [2, 3])
        @test_throws ArgumentError partialtrace(rand(6, 6), 1, [2, 4])
        @test_throws ArgumentError partialtrace(rand(3, 4), 1, [2, 3])
    end
end

@add_problem sdp function sdp_Real_Variables_with_complex_equality_constraints(
    handle_problem!,
    ::Val{test},
    atol,
    rtol,
    ::Type{T},
) where {T,test}
    n = 10 # variable dimension (parameter)
    m = 5 # number of constraints (parameter)
    xo = rand(n)
    A = randn(m, n) + im * randn(m, n)
    b = A * xo
    x = Variable(n)
    p1 = minimize(sum(x), A * x == b, x >= 0; numeric_type = T)

    handle_problem!(p1)

    if test
        x1 = copy(evaluate(x))
    end

    p2 = minimize(
        sum(x),
        real(A) * x == real(b),
        imag(A) * x == imag(b),
        x >= 0;
        numeric_type = T,
    )

    handle_problem!(p2)
    if test
        x2 = evaluate(x)
        @test x1 == x2
    end
end

@add_problem sdp function sdp_Complex_Variable_with_complex_equality_constraints(
    handle_problem!,
    ::Val{test},
    atol,
    rtol,
    ::Type{T},
) where {T,test}
    n = 10 # variable dimension (parameter)
    m = 5 # number of constraints (parameter)
    xo = rand(n) + im * rand(n)
    A = randn(m, n) + im * randn(m, n)
    b = A * xo
    x = ComplexVariable(n)
    p1 = minimize(
        real(sum(x)),
        A * x == b,
        real(x) >= 0,
        imag(x) >= 0;
        numeric_type = T,
    )

    handle_problem!(p1)

    if test
        x1 = evaluate(x)
    end

    xr = Variable(n)
    xi = Variable(n)
    p2 = minimize(
        sum(xr),
        real(A) * xr - imag(A) * xi == real(b),
        imag(A) * xr + real(A) * xi == imag(b),
        xr >= 0,
        xi >= 0;
        numeric_type = T,
    )

    handle_problem!(p2)

    if test
        @test real(x1) ≈ evaluate(xr) atol = atol rtol = rtol
        @test imag(x1) ≈ evaluate(xi) atol = atol rtol = rtol
    end
end

@add_problem sdp function sdp_Issue_198(
    handle_problem!,
    ::Val{test},
    atol,
    rtol,
    ::Type{T},
) where {T,test}
    ρ = HermitianSemidefinite(2)
    constraints = [ρ == [1.0 0.0; 0.0 1.0]]
    p = satisfy(constraints; numeric_type = T)

    handle_problem!(p)
    if test
        @test p.status == MOI.OPTIMAL
        @test evaluate(ρ) ≈ [1.0 0.0; 0.0 1.0] atol = atol rtol = rtol
        @test p.optval === nothing
    end
end

@add_problem sdp function sdp_socp_norm2_atom(
    handle_problem!,
    ::Val{test},
    atol,
    rtol,
    ::Type{T},
) where {T,test}
    a = 2 + 4im
    x = ComplexVariable()
    objective = norm2(a + (-x))
    c1 = real(x) >= 0
    p = minimize(objective, c1; numeric_type = T)

    handle_problem!(p)
    if test
        @test p.optval ≈ 0 atol = atol rtol = rtol
        @test evaluate(objective) ≈ 0 atol = atol rtol = rtol

        @test real(evaluate(x)) ≈ real(a) atol = atol rtol = rtol
        @test imag(evaluate(x)) ≈ imag(a) atol = atol rtol = rtol
    end
end

@add_problem sdp function sdp_socp_sumsquares_atom(
    handle_problem!,
    ::Val{test},
    atol,
    rtol,
    ::Type{T},
) where {T,test}
    a = [2 + 4im; 4 + 6im]
    x = ComplexVariable(2)
    objective = sumsquares(a - x)
    c1 = real(x) >= 0
    p = minimize(objective, c1; numeric_type = T)

    handle_problem!(p)
    if test
        @test p.optval ≈ 0 atol = atol rtol = rtol
        @test evaluate(objective) ≈ 0.0 atol = atol rtol = rtol

        @test real.(evaluate(x)) ≈ real.(a) atol = atol rtol = rtol
        @test imag.(evaluate(x)) ≈ imag.(a) atol = atol rtol = rtol
    end
end

@add_problem sdp function sdp_socp_abs_atom(
    handle_problem!,
    ::Val{test},
    atol,
    rtol,
    ::Type{T},
) where {T,test}
    a = 5 - 4im
    x = ComplexVariable()
    objective = abs(a - x)
    c1 = real(x) >= 0
    p = minimize(objective, c1; numeric_type = T)

    handle_problem!(p)
    if test
        @test p.optval ≈ 0 atol = atol rtol = rtol
        @test evaluate(objective) ≈ 0.0 atol = atol rtol = rtol

        @test real(evaluate(x)) ≈ real(a) atol = atol rtol = rtol
        @test imag(evaluate(x)) ≈ imag(a) atol = atol rtol = rtol
    end
end

@add_problem sdp function sdp_Complex_Semidefinite_constraint(
    handle_problem!,
    ::Val{test},
    atol,
    rtol,
    ::Type{T},
) where {T,test}
    n = 10
    A = rand(n, n) + im * rand(n, n)
    A = A + A' # now A is hermitian
    x = ComplexVariable(n, n)
    objective = sumsquares(A - x)
    c1 = isposdef(x)
    p = minimize(objective, c1; numeric_type = T)

    handle_problem!(p)
    # test that X is approximately equal to posA:
    l, v = LinearAlgebra.eigen(A)
    posA = v * Diagonal(max.(l, 0)) * v'

    if test
        @test real.(evaluate(x)) ≈ real.(posA) atol = atol rtol = rtol
        @test imag.(evaluate(x)) ≈ imag.(posA) atol = atol rtol = rtol
    end
end

@add_problem sdp function sdp_relative_entropy(
    handle_problem!,
    ::Val{test},
    atol,
    rtol,
    ::Type{T},
) where {T,test}
    for cplx in [false, true]
        n = 4
        X = eye(n)
        if cplx
            Y = ComplexVariable(n, n)
        else
            Y = Variable(n, n)
        end

        c1 = Convex.Constraint((eye(n), X, Y), RelativeEntropyEpiConeSquare(n))
        objective = real(tr(Y))
        p = minimize(objective, c1; numeric_type = T)

        handle_problem!(p)

        if test
            @test evaluate(Y) ≈ eye(n) * exp(-1) atol = atol rtol = rtol
        end
    end
end

@add_problem sdp function sdp_geom_mean_hypocone_real_0(
    handle_problem!,
    ::Val{test},
    atol,
    rtol,
    ::Type{T},
) where {T,test}
    n = 4
    A = Variable(n, n)
    B = randn(n, n)
    B = B * B' # now A is positive semidefinite
    B += 0.2 * LinearAlgebra.I # prevent numerical instability

    c1 = Convex.Constraint(
        (eye(n), A, B),
        GeometricMeanHypoConeSquare(0 // 1, n),
    )
    objective = tr(A)
    p = minimize(objective, c1; numeric_type = T)

    handle_problem!(p)

    if test
        @test evaluate(A) ≈ eye(n) atol = atol rtol = rtol
    end
end

@add_problem sdp function sdp_geom_mean_hypocone_real_1(
    handle_problem!,
    ::Val{test},
    atol,
    rtol,
    ::Type{T},
) where {T,test}
    n = 4
    A = Variable(n, n)
    B = randn(n, n)
    B = B * B' # now A is positive semidefinite
    B += 0.2 * LinearAlgebra.I # prevent numerical instability

    c1 = Convex.Constraint(
        (eye(n), B, A),
        GeometricMeanHypoConeSquare(1 // 1, n),
    )
    objective = tr(A)
    p = minimize(objective, c1; numeric_type = T)

    handle_problem!(p)

    if test
        @test evaluate(A) ≈ eye(n) atol = atol rtol = rtol
    end
end

@add_problem sdp function sdp_geom_mean_hypocone_real_1_2(
    handle_problem!,
    ::Val{test},
    atol,
    rtol,
    ::Type{T},
) where {T,test}
    n = 4
    A = randn(n, n)
    A = A * A' # now A is positive semidefinite
    A += 0.2 * LinearAlgebra.I # prevent numerical instability
    B = Variable(n, n)

    c1 = Convex.Constraint(
        (eye(n), A, B),
        GeometricMeanHypoConeSquare(1 // 2, n),
    )
    objective = tr(B)
    p = minimize(objective, c1; numeric_type = T)

    handle_problem!(p)

    if test
        @test evaluate(B) ≈ A^-1 atol = atol rtol = rtol
    end
end

@add_problem sdp function sdp_geom_mean_hypocone_real_3_8(
    handle_problem!,
    ::Val{test},
    atol,
    rtol,
    ::Type{T},
) where {T,test}
    n = 4
    A = randn(n, n)
    A = A * A' # now A is positive semidefinite
    A += 0.2 * LinearAlgebra.I # prevent numerical instability
    B = Variable(n, n)

    c1 = Convex.Constraint(
        (eye(n), A, B),
        GeometricMeanHypoConeSquare(3 // 8, n),
    )
    objective = tr(B)
    p = minimize(objective, c1; numeric_type = T)

    handle_problem!(p)

    # A #_t B = I  =>  B = A^(1-1/t)
    if test
        @test evaluate(B) ≈ A^(-5 // 3) atol = atol rtol = rtol
    end
end

@add_problem sdp function sdp_geom_mean_hypocone_real_3_5(
    handle_problem!,
    ::Val{test},
    atol,
    rtol,
    ::Type{T},
) where {T,test}
    n = 4
    A = randn(n, n)
    A = A * A' # now A is positive semidefinite
    A += 0.2 * LinearAlgebra.I # prevent numerical instability
    B = Variable(n, n)

    c1 = Convex.Constraint(
        (eye(n), A, B),
        GeometricMeanHypoConeSquare(3 // 5, n),
    )
    objective = tr(B)
    p = minimize(objective, c1; numeric_type = T)

    handle_problem!(p)

    # A #_t B = I  =>  B = A^(1-1/t)
    if test
        @test evaluate(B) ≈ A^(-2 // 3) atol = atol rtol = rtol
    end
end

@add_problem sdp function sdp_geom_mean_hypocone_cplx_1_2(
    handle_problem!,
    ::Val{test},
    atol,
    rtol,
    ::Type{T},
) where {T,test}
    n = 4
    A = randn(ComplexF64, n, n)
    A = A * A' # now A is positive semidefinite
    A += 0.2 * LinearAlgebra.I # prevent numerical instability
    B = ComplexVariable(n, n)

    c1 = Convex.Constraint(
        (eye(n), A, B),
        GeometricMeanHypoConeSquare(1 // 2, n),
    )
    objective = real(tr(B))
    p = minimize(objective, c1; numeric_type = T)

    handle_problem!(p)

    if test
        @test real.(evaluate(B)) ≈ real.(A^-1) atol = atol rtol = rtol
        @test imag.(evaluate(B)) ≈ imag.(A^-1) atol = atol rtol = rtol
    end
end

@add_problem sdp function sdp_geom_mean_hypocone_cplx_3_8(
    handle_problem!,
    ::Val{test},
    atol,
    rtol,
    ::Type{T},
) where {T,test}
    n = 4
    A = randn(ComplexF64, n, n)
    A = A * A' # now A is positive semidefinite
    A += 0.2 * LinearAlgebra.I # prevent numerical instability
    B = ComplexVariable(n, n)

    c1 = Convex.Constraint(
        (eye(n), A, B),
        GeometricMeanHypoConeSquare(3 // 8, n),
    )
    objective = real(tr(B))
    p = minimize(objective, c1; numeric_type = T)

    handle_problem!(p)

    # A #_t B = I  =>  B = A^(1-1/t)
    if test
        @test real.(evaluate(B)) ≈ real.(A^(-5 // 3)) atol = atol rtol = rtol
        @test imag.(evaluate(B)) ≈ imag.(A^(-5 // 3)) atol = atol rtol = rtol
    end
end

@add_problem sdp function sdp_geom_mean_hypocone_cplx_3_5(
    handle_problem!,
    ::Val{test},
    atol,
    rtol,
    ::Type{T},
) where {T,test}
    n = 4
    A = randn(ComplexF64, n, n)
    A = A * A' # now A is positive semidefinite
    A += 0.2 * LinearAlgebra.I # prevent numerical instability
    B = ComplexVariable(n, n)

    c1 = Convex.Constraint(
        (eye(n), A, B),
        GeometricMeanHypoConeSquare(3 // 5, n),
    )
    objective = real(tr(B))
    p = minimize(objective, c1; numeric_type = T)

    handle_problem!(p)

    # A #_t B = I  =>  B = A^(1-1/t)
    if test
        @test real.(evaluate(B)) ≈ real.(A^(-2 // 3)) atol = atol rtol = rtol
        @test imag.(evaluate(B)) ≈ imag.(A^(-2 // 3)) atol = atol rtol = rtol
    end
end

@add_problem sdp function sdp_geom_mean_hypocone_fullhyp(
    handle_problem!,
    ::Val{test},
    atol,
    rtol,
    ::Type{T},
) where {T,test}
    # These matrices satisfy A^(1/2) ⪰ B and so (A,eye(2),B) \in hyp_{1/2}
    # A and B however do not satisfy A ⪰ B^2 and so [A B; B eye(2)] not
    # psd, and using fullhyp=false will give an infeasible SDP.
    A = [6.25 0; 0 16]
    B = Variable(2, 2)
    S = Semidefinite(2)

    p = minimize(
        0,
        [
            B == [2 1; 1 2],
            Convex.Constraint(
                (B, A, eye(2)),
                GeometricMeanHypoConeSquare(1 // 2, 2, false),
            ),
        ];
        numeric_type = T,
    )
    handle_problem!(p)
    if test
        @test p.status == MOI.INFEASIBLE
    end

    p = minimize(
        0,
        [
            B == [2 1; 1 2],
            Convex.Constraint(
                (B, A, eye(2)),
                GeometricMeanHypoConeSquare(1 // 2, 2),
            ),
        ];
        numeric_type = T,
    )
    handle_problem!(p)
    if test
        @test p.status in (MOI.OPTIMAL, MOI.ALMOST_OPTIMAL)
    end
end

@add_problem sdp function sdp_geom_mean_epicone_real_neg1_optA(
    handle_problem!,
    ::Val{test},
    atol,
    rtol,
    ::Type{T},
) where {T,test}
    n = 3
    B = randn(n, n)
    B = B * B' # now B is positive semidefinite
    B += 0.2 * LinearAlgebra.I # prevent numerical instability
    A = Variable(n, n)

    c1 = Convex.Constraint(
        (eye(n), A, B),
        GeometricMeanEpiConeSquare(-1 // 1, n),
    )
    objective = tr(A)
    p = maximize(objective, c1; numeric_type = T)

    handle_problem!(p)

    # A #_t B = I  =>  B = A^(1-1/t)
    if test
        @test B ≈ evaluate(A)^2 atol = atol rtol = rtol
    end
end

@add_problem sdp function sdp_geom_mean_epicone_real_neg1_optB(
    handle_problem!,
    ::Val{test},
    atol,
    rtol,
    ::Type{T},
) where {T,test}
    n = 3
    A = randn(n, n)
    A = A * A' # now A is positive semidefinite
    A += 0.2 * LinearAlgebra.I # prevent numerical instability
    A /= tr(A) # solver has problems if B is large
    B = Variable(n, n)

    c1 = Convex.Constraint(
        (eye(n), A, B),
        GeometricMeanEpiConeSquare(-1 // 1, n),
    )
    objective = tr(B)
    p = minimize(objective, c1; numeric_type = T)

    handle_problem!(p)

    # A #_t B = I  =>  B = A^(1-1/t)
    if test
        @test evaluate(B) ≈ A^2 atol = atol rtol = rtol
    end
end

@add_problem sdp function sdp_geom_mean_epicone_real_neg3_5_optA(
    handle_problem!,
    ::Val{test},
    atol,
    rtol,
    ::Type{T},
) where {T,test}
    n = 3
    B = randn(n, n)
    B = B * B' # now B is positive semidefinite
    B += 0.2 * LinearAlgebra.I # prevent numerical instability
    A = Variable(n, n)

    c1 = Convex.Constraint(
        (eye(n), A, B),
        GeometricMeanEpiConeSquare(-3 // 5, n),
    )
    objective = tr(A)
    p = maximize(objective, c1; numeric_type = T)

    handle_problem!(p)

    # A #_t B = I  =>  B = A^(1-1/t)
    if test
        @test B ≈ evaluate(A)^(8 // 3) atol = atol rtol = rtol
    end
end

@add_problem sdp function sdp_geom_mean_epicone_real_neg3_5_optB(
    handle_problem!,
    ::Val{test},
    atol,
    rtol,
    ::Type{T},
) where {T,test}
    n = 3
    A = randn(n, n)
    A = A * A' # now A is positive semidefinite
    A += 0.2 * LinearAlgebra.I # prevent numerical instability
    A /= tr(A) # solver has problems if B is large
    B = Variable(n, n)

    c1 = Convex.Constraint(
        (eye(n), A, B),
        GeometricMeanEpiConeSquare(-3 // 5, n),
    )
    objective = tr(B)
    p = minimize(objective, c1; numeric_type = T)

    handle_problem!(p)

    # A #_t B = I  =>  B = A^(1-1/t)
    if test
        @test evaluate(B) ≈ A^(8 // 3) atol = atol rtol = rtol
    end
end

@add_problem sdp function sdp_geom_mean_epicone_real_8_5_optA(
    handle_problem!,
    ::Val{test},
    atol,
    rtol,
    ::Type{T},
) where {T,test}
    n = 3
    B = randn(n, n)
    B = B * B' # now B is positive semidefinite
    B += 0.2 * LinearAlgebra.I # prevent numerical instability
    B /= tr(B) # solver has problems if B is large
    A = Variable(n, n)

    c1 =
        Convex.Constraint((eye(n), A, B), GeometricMeanEpiConeSquare(8 // 5, n))
    objective = tr(A)
    p = minimize(objective, c1; numeric_type = T)

    handle_problem!(p)

    # A #_t B = I  =>  B = A^(1-1/t)
    if test
        @test B ≈ evaluate(A)^(3 // 8) atol = atol rtol = rtol
    end
end

@add_problem sdp function sdp_geom_mean_epicone_real_8_5_optB(
    handle_problem!,
    ::Val{test},
    atol,
    rtol,
    ::Type{T},
) where {T,test}
    n = 3
    A = randn(n, n)
    A = A * A' # now A is positive semidefinite
    A += 0.2 * LinearAlgebra.I # prevent numerical instability
    B = Variable(n, n)

    c1 =
        Convex.Constraint((eye(n), A, B), GeometricMeanEpiConeSquare(8 // 5, n))
    objective = tr(B)
    p = maximize(objective, c1; numeric_type = T)

    handle_problem!(p)

    # A #_t B = I  =>  B = A^(1-1/t)
    if test
        @test evaluate(B) ≈ A^(3 // 8) atol = atol rtol = rtol
    end
end

@add_problem sdp function sdp_quantum_relative_entropy_argcheck(
    handle_problem!,
    ::Val{test},
    atol,
    rtol,
    ::Type{T},
) where {T,test}
    if test
        @test_throws DimensionMismatch quantum_relative_entropy(
            zeros(2, 3),
            zeros(2, 3),
        )
        @test_throws DimensionMismatch quantum_relative_entropy(
            zeros(2, 3),
            Variable(2, 3),
        )
        @test_throws DimensionMismatch quantum_relative_entropy(
            Variable(2, 3),
            zeros(2, 3),
        )
        @test_throws DimensionMismatch quantum_relative_entropy(
            Variable(2, 3),
            Variable(2, 3),
        )

        @test_throws DimensionMismatch quantum_relative_entropy(
            zeros(2, 2),
            zeros(3, 3),
        )
        @test_throws DimensionMismatch quantum_relative_entropy(
            zeros(2, 2),
            Variable(3, 3),
        )
        @test_throws DimensionMismatch quantum_relative_entropy(
            Variable(2, 2),
            zeros(3, 3),
        )
        @test_throws DimensionMismatch quantum_relative_entropy(
            Variable(2, 2),
            Variable(2, 3),
        )

        z = zeros(2, 2)
        v = Variable(2, 2)
        nh = [1 1; 0 1] # not hermitian
        np = [1 2; 2 1] # not positive semidefinite
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

function sdp_quantum_relative_entropy_impl(
    handle_problem!,
    ::Val{test},
    atol,
    rtol,
    ::Type{T},
    lowrank::Bool,
    mode::Integer,
) where {T,test}
    # Too slow for n>3.  Anything that can be done to speed it up?
    n = 3
    X = randn(ComplexF64, n, n)
    if lowrank
        X[:, 1] .= 0
    end
    X = X * X' # now X is positive semidefinite
    X /= tr(X) # make it a quantum state

    A = ComplexVariable(n, n)
    B = ComplexVariable(n, n)
    constraints = Array{Constraint,1}([A == X])
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
        @test real.(evaluate(A)) ≈ real.(X) atol = atol rtol = rtol
        @test imag.(evaluate(A)) ≈ imag.(X) atol = atol rtol = rtol
        if mode != 5
            @test real.(evaluate(B)) ≈ real.(X) atol = atol rtol = rtol
            @test imag.(evaluate(B)) ≈ imag.(X) atol = atol rtol = rtol
            @test p.optval ≈ 0 atol = atol rtol = rtol
        end
        if mode == 1
            @test p.optval ≈ evaluate(quantum_relative_entropy(A, B)) atol =
                atol rtol = rtol
        elseif mode == 2
            @test p.optval ≈ evaluate(quantum_relative_entropy(X, B)) atol =
                atol rtol = rtol
        elseif mode == 3
            @test p.optval ≈ evaluate(quantum_relative_entropy(B, A)) atol =
                atol rtol = rtol
        elseif mode == 4
            @test p.optval ≈ evaluate(quantum_relative_entropy(B, X)) atol =
                atol rtol = rtol
        elseif mode == 5
            # Satisfiability problem
            @test p.optval ≈ 0 atol = atol rtol = rtol
        end
    end
end

@add_problem sdp function sdp_quantum_relative_entropy1_fullrank(
    handle_problem!,
    ::Val{test},
    atol,
    rtol,
    ::Type{T},
) where {T,test}
    return sdp_quantum_relative_entropy_impl(
        handle_problem!,
        Val(test),
        atol,
        rtol,
        T,
        false,
        1,
    )
end

@add_problem sdp function sdp_quantum_relative_entropy2_fullrank(
    handle_problem!,
    ::Val{test},
    atol,
    rtol,
    ::Type{T},
) where {T,test}
    return sdp_quantum_relative_entropy_impl(
        handle_problem!,
        Val(test),
        atol,
        rtol,
        T,
        false,
        2,
    )
end

@add_problem sdp function sdp_quantum_relative_entropy3_fullrank(
    handle_problem!,
    ::Val{test},
    atol,
    rtol,
    ::Type{T},
) where {T,test}
    return sdp_quantum_relative_entropy_impl(
        handle_problem!,
        Val(test),
        atol,
        rtol,
        T,
        false,
        3,
    )
end

@add_problem sdp function sdp_quantum_relative_entropy4_fullrank(
    handle_problem!,
    ::Val{test},
    atol,
    rtol,
    ::Type{T},
) where {T,test}
    return sdp_quantum_relative_entropy_impl(
        handle_problem!,
        Val(test),
        atol,
        rtol,
        T,
        false,
        4,
    )
end

@add_problem sdp function sdp_quantum_relative_entropy5_fullrank(
    handle_problem!,
    ::Val{test},
    atol,
    rtol,
    ::Type{T},
) where {T,test}
    return sdp_quantum_relative_entropy_impl(
        handle_problem!,
        Val(test),
        atol,
        rtol,
        T,
        false,
        5,
    )
end

@add_problem sdp function sdp_quantum_relative_entropy1_lowrank(
    handle_problem!,
    ::Val{test},
    atol,
    rtol,
    ::Type{T},
) where {T,test}
    return sdp_quantum_relative_entropy_impl(
        handle_problem!,
        Val(test),
        atol,
        rtol,
        T,
        true,
        1,
    )
end

@add_problem sdp function sdp_quantum_relative_entropy2_lowrank(
    handle_problem!,
    ::Val{test},
    atol,
    rtol,
    ::Type{T},
) where {T,test}
    return sdp_quantum_relative_entropy_impl(
        handle_problem!,
        Val(test),
        atol,
        rtol,
        T,
        true,
        2,
    )
end

@add_problem sdp function sdp_quantum_relative_entropy3_lowrank(
    handle_problem!,
    ::Val{test},
    atol,
    rtol,
    ::Type{T},
) where {T,test}
    return sdp_quantum_relative_entropy_impl(
        handle_problem!,
        Val(test),
        atol,
        rtol,
        T,
        true,
        3,
    )
end

@add_problem sdp function sdp_quantum_relative_entropy4_lowrank(
    handle_problem!,
    ::Val{test},
    atol,
    rtol,
    ::Type{T},
) where {T,test}
    return sdp_quantum_relative_entropy_impl(
        handle_problem!,
        Val(test),
        atol,
        rtol,
        T,
        true,
        4,
    )
end

@add_problem sdp function sdp_quantum_relative_entropy5_lowrank(
    handle_problem!,
    ::Val{test},
    atol,
    rtol,
    ::Type{T},
) where {T,test}
    return sdp_quantum_relative_entropy_impl(
        handle_problem!,
        Val(test),
        atol,
        rtol,
        T,
        true,
        5,
    )
end

@add_problem sdp function sdp_quantum_relative_entropy_const(
    handle_problem!,
    ::Val{test},
    atol,
    rtol,
    ::Type{T},
) where {T,test}
    if test
        @test quantum_relative_entropy(
            diagm(0 => [0.5, 0.5]),
            diagm(0 => [0.5, 0.5]),
        ) ≈ 0 atol = atol rtol = rtol
        @test quantum_relative_entropy(
            diagm(0 => [1.0, 0.0]),
            diagm(0 => [1.0, 0.0]),
        ) ≈ 0 atol = atol rtol = rtol
        @test quantum_relative_entropy(
            diagm(0 => [1.0, 0.0]),
            diagm(0 => [0.5, 0.5]),
        ) ≈ log(2) atol = atol rtol = rtol
        @test isinf(
            quantum_relative_entropy(
                diagm(0 => [0.5, 0.5]),
                diagm(0 => [1.0, 0.0]),
            ),
        )
    end
end

@add_problem sdp function sdp_quantum_channel_capacity(
    handle_problem!,
    ::Val{test},
    atol,
    rtol,
    ::Type{T},
) where {T,test}
    # Capacity of classical-quantum channel.  See chapter 3 of arxiv:1705.06671,
    # Fawzi, Fawzi "Efficient optimization of the quantum relative entropy".
    # For channels with two symbols, both pure states, we have an exact formula to compare against.

    n = 3

    for cplx in [false, true]
        if cplx
            ψ1 = LinearAlgebra.normalize(randn(ComplexF64, n))
            ψ2 = LinearAlgebra.normalize(randn(ComplexF64, n))
        else
            ψ1 = LinearAlgebra.normalize(randn(Float64, n))
            ψ2 = LinearAlgebra.normalize(randn(Float64, n))
        end
        ρ1 = ψ1 * ψ1'
        ρ2 = ψ2 * ψ2'

        p1 = Variable()
        p2 = Variable()

        objective =
            (
                quantum_entropy(p1 * ρ1 + p2 * ρ2) - p1 * quantum_entropy(ρ1) -
                p2 * quantum_entropy(ρ2)
            ) / log(2)
        constraints = [p1 >= 0, p2 >= 0, p1 + p2 == 1]
        p = maximize(objective, constraints; numeric_type = T)
        handle_problem!(p)

        ϵ = abs(ψ1' * ψ2)
        q = (1 + ϵ) / 2
        r = [q, 1 - q]
        v = -r' * log2.(r)

        if test
            @test p.optval ≈ v atol = atol rtol = rtol
            @test evaluate(objective) ≈ v atol = atol rtol = rtol

            @test_throws DimensionMismatch quantum_entropy(Variable(2, 3))
        end
    end
end

@add_problem sdp function sdp_trace_mpower_argcheck(
    handle_problem!,
    ::Val{test},
    atol,
    rtol,
    ::Type{T},
) where {T,test}
    if test
        @test_throws DimensionMismatch trace_mpower(
            zeros(2, 3),
            1 // 2,
            zeros(2, 3),
        )
        @test_throws DimensionMismatch trace_mpower(
            Variable(2, 3),
            1 // 2,
            zeros(2, 3),
        )

        @test_throws DimensionMismatch trace_mpower(
            zeros(2, 2),
            1 // 2,
            zeros(3, 3),
        )
        @test_throws DimensionMismatch trace_mpower(
            Variable(2, 2),
            1 // 2,
            zeros(3, 3),
        )

        @test_throws DomainError trace_mpower(
            Variable(3, 3),
            5 // 2,
            zeros(3, 3),
        )

        nh = [1 1; 0 1] # not hermitian
        np = [1 2; 2 1] # not positive semidefinite
        @test_throws DomainError trace_mpower(Variable(2, 2), 1 // 2, nh)
        @test_throws DomainError trace_mpower(Variable(2, 2), 1 // 2, np)
    end
end

@add_problem sdp function sdp_trace_mpower_real_2_3(
    handle_problem!,
    ::Val{test},
    atol,
    rtol,
    ::Type{T},
) where {T,test}
    # Too slow for n>3.  Anything that can be done to speed it up?
    n = 3
    A = randn(n, n)
    A = A * A' # now A is positive semidefinite
    C = randn(n, n)
    C = C * C' # now C is positive semidefinite
    t = 2 // 3

    B = Semidefinite(n, n)

    constraints = [B ⪯ A]
    objective = trace_mpower(B, t, C)
    p = maximize(objective, constraints; numeric_type = T)

    handle_problem!(p)

    if test
        @test real.(evaluate(B)) ≈ real.(A) atol = atol rtol = rtol
        @test imag.(evaluate(B)) ≈ imag.(A) atol = atol rtol = rtol
        @test p.optval ≈ tr(C * A^t) atol = atol rtol = rtol
        @test p.optval ≈ trace_mpower(A, t, C) atol = atol rtol = rtol
        @test p.optval ≈ evaluate(trace_mpower(B, t, C)) atol = atol rtol = rtol
    end
end

@add_problem sdp function sdp_trace_mpower_cplx_2_3(
    handle_problem!,
    ::Val{test},
    atol,
    rtol,
    ::Type{T},
) where {T,test}
    # Too slow for n>3.  Anything that can be done to speed it up?
    n = 3
    A = randn(ComplexF64, n, n)
    A = A * A' # now A is positive semidefinite
    C = randn(ComplexF64, n, n)
    C = C * C' # now C is positive semidefinite
    t = 2 // 3

    B = HermitianSemidefinite(n, n)

    constraints = [B ⪯ A]
    objective = trace_mpower(B, t, C)
    p = maximize(objective, constraints; numeric_type = T)

    handle_problem!(p)

    if test
        @test real.(evaluate(B)) ≈ real.(A) atol = atol rtol = rtol
        @test imag.(evaluate(B)) ≈ imag.(A) atol = atol rtol = rtol
        @test p.optval ≈ tr(C * A^t) atol = atol rtol = rtol
        @test p.optval ≈ trace_mpower(A, t, C) atol = atol rtol = rtol
        @test p.optval ≈ evaluate(trace_mpower(B, t, C)) atol = atol rtol = rtol
    end
end

@add_problem sdp function sdp_trace_mpower_real_5_4(
    handle_problem!,
    ::Val{test},
    atol,
    rtol,
    ::Type{T},
) where {T,test}
    # Too slow for n>3.  Anything that can be done to speed it up?
    n = 3
    A = randn(n, n)
    A = A * A' # now A is positive semidefinite
    C = randn(n, n)
    C = C * C' # now C is positive semidefinite
    t = 5 // 4

    B = Semidefinite(n, n)

    constraints = [
        # note: B ⪰ A doesn't work for this unit test because trace_mpower is not operator monotone
        # unless t ∈ [0,1].
        B == A,
    ]
    objective = trace_mpower(B, t, C)
    p = minimize(objective, constraints; numeric_type = T)

    handle_problem!(p)

    if test
        @test real.(evaluate(B)) ≈ real.(A) atol = atol rtol = rtol
        @test imag.(evaluate(B)) ≈ imag.(A) atol = atol rtol = rtol
        @test p.optval ≈ tr(C * A^t) atol = atol rtol = rtol
        @test p.optval ≈ trace_mpower(A, t, C) atol = atol rtol = rtol
        @test p.optval ≈ evaluate(trace_mpower(B, t, C)) atol = atol rtol = rtol
    end
end

@add_problem sdp function sdp_trace_mpower_cplx_5_4(
    handle_problem!,
    ::Val{test},
    atol,
    rtol,
    ::Type{T},
) where {T,test}
    # Too slow for n>3.  Anything that can be done to speed it up?
    n = 3
    A = randn(ComplexF64, n, n)
    A = A * A' # now A is positive semidefinite
    C = randn(ComplexF64, n, n)
    C = C * C' # now C is positive semidefinite
    t = 5 // 4

    B = HermitianSemidefinite(n, n)

    constraints = [
        # note: B ⪰ A doesn't work for this unit test because trace_mpower is not operator monotone
        # unless t ∈ [0,1].
        B == A,
    ]
    objective = trace_mpower(B, t, C)
    p = minimize(objective, constraints; numeric_type = T)

    handle_problem!(p)

    if test
        @test real.(evaluate(B)) ≈ real.(A) atol = atol rtol = rtol
        @test imag.(evaluate(B)) ≈ imag.(A) atol = atol rtol = rtol
        @test p.optval ≈ tr(C * A^t) atol = atol rtol = rtol
        @test p.optval ≈ trace_mpower(A, t, C) atol = atol rtol = rtol
        @test p.optval ≈ evaluate(trace_mpower(B, t, C)) atol = atol rtol = rtol
    end
end

@add_problem sdp function sdp_trace_mpower_real_neg1_4(
    handle_problem!,
    ::Val{test},
    atol,
    rtol,
    ::Type{T},
) where {T,test}
    # Too slow for n>3.  Anything that can be done to speed it up?
    n = 3
    A = randn(n, n)
    A = A * A' # now A is positive semidefinite
    A += 0.2 * LinearAlgebra.I # prevent numerical instability
    C = randn(n, n)
    C = C * C' # now C is positive semidefinite
    t = -1 // 4

    B = Semidefinite(n, n)

    constraints = [
        # -trace_mpower is operator monotone for t ∈ [-1, 0].
        B ⪯ A,
    ]
    objective = trace_mpower(B, t, C)
    p = minimize(objective, constraints; numeric_type = T)

    handle_problem!(p)

    if test
        @test real.(evaluate(B)) ≈ real.(A) atol = atol rtol = rtol
        @test imag.(evaluate(B)) ≈ imag.(A) atol = atol rtol = rtol
        @test p.optval ≈ tr(C * A^t) atol = atol rtol = rtol
        @test p.optval ≈ trace_mpower(A, t, C) atol = atol rtol = rtol
        @test p.optval ≈ evaluate(trace_mpower(B, t, C)) atol = atol rtol = rtol
    end
end

@add_problem sdp function sdp_trace_mpower_cplx_neg1_4(
    handle_problem!,
    ::Val{test},
    atol,
    rtol,
    ::Type{T},
) where {T,test}
    # Too slow for n>3.  Anything that can be done to speed it up?
    n = 3
    A = randn(ComplexF64, n, n)
    A = A * A' # now A is positive semidefinite
    A += 0.2 * LinearAlgebra.I # prevent numerical instability
    C = randn(ComplexF64, n, n)
    C = C * C' # now C is positive semidefinite
    t = -1 // 4

    B = HermitianSemidefinite(n, n)

    constraints = [
        # -trace_mpower is operator monotone for t ∈ [-1, 0].
        B ⪯ A,
    ]
    objective = trace_mpower(B, t, C)
    p = minimize(objective, constraints; numeric_type = T)

    handle_problem!(p)

    if test
        @test real.(evaluate(B)) ≈ real.(A) atol = atol rtol = rtol
        @test imag.(evaluate(B)) ≈ imag.(A) atol = atol rtol = rtol
        @test p.optval ≈ tr(C * A^t) atol = atol rtol = rtol
        @test p.optval ≈ trace_mpower(A, t, C) atol = atol rtol = rtol
        @test p.optval ≈ evaluate(trace_mpower(B, t, C)) atol = atol rtol = rtol
    end
end

@add_problem sdp function sdp_trace_logm_argcheck(
    handle_problem!,
    ::Val{test},
    atol,
    rtol,
    ::Type{T},
) where {T,test}
    if test
        @test_throws DimensionMismatch trace_logm(zeros(2, 3), zeros(2, 3))
        @test_throws DimensionMismatch trace_logm(Variable(2, 3), zeros(2, 3))

        @test_throws DimensionMismatch trace_logm(zeros(2, 2), zeros(3, 3))
        @test_throws DimensionMismatch trace_logm(Variable(2, 2), zeros(3, 3))

        nh = [1 1; 0 1] # not hermitian
        np = [1 2; 2 1] # not positive semidefinite
        @test_throws DomainError trace_logm(Variable(2, 2), nh)
        @test_throws DomainError trace_logm(Variable(2, 2), np)
    end
end

@add_problem sdp function sdp_trace_logm_real(
    handle_problem!,
    ::Val{test},
    atol,
    rtol,
    ::Type{T},
) where {T,test}
    n = 3
    C = randn(Float64, n, n)
    C = C * C' # now C is positive semidefinite

    X = Variable(n, n)

    constraints = [X ⪯ eye(n)]
    objective = trace_logm(X, C)
    p = maximize(objective, constraints; numeric_type = T)

    handle_problem!(p)

    if test
        @test real.(evaluate(X)) ≈ real.(eye(n)) atol = atol rtol = rtol
        @test imag.(evaluate(X)) ≈ imag.(eye(n)) atol = atol rtol = rtol
        @test p.optval ≈ tr(C * log(evaluate(X))) atol = atol rtol = rtol
        @test p.optval ≈ trace_logm(evaluate(X), C) atol = atol rtol = rtol
        @test p.optval ≈ evaluate(trace_logm(X, C)) atol = atol rtol = rtol
    end
end

@add_problem sdp function sdp_trace_logm_cplx(
    handle_problem!,
    ::Val{test},
    atol,
    rtol,
    ::Type{T},
) where {T,test}
    n = 3
    C = randn(ComplexF64, n, n)
    C = C * C' # now C is positive semidefinite

    X = ComplexVariable(n, n)

    constraints = [X ⪯ eye(n)]
    objective = trace_logm(X, C)
    p = maximize(objective, constraints; numeric_type = T)

    handle_problem!(p)

    if test
        @test real.(evaluate(X)) ≈ real.(eye(n)) atol = atol rtol = rtol
        @test imag.(evaluate(X)) ≈ imag.(eye(n)) atol = atol rtol = rtol
        @test p.optval ≈ tr(C * log(evaluate(X))) atol = atol rtol = rtol
        @test p.optval ≈ trace_logm(evaluate(X), C) atol = atol rtol = rtol
        @test p.optval ≈ evaluate(trace_logm(X, C)) atol = atol rtol = rtol
    end
end

@add_problem sdp function sdp_lieb_ando(
    handle_problem!,
    ::Val{test},
    atol,
    rtol,
    ::Type{T},
) where {T,test}
    if test
        @test_throws DomainError lieb_ando(
            zeros(3, 3),
            zeros(3, 3),
            eye(3),
            -3 // 2,
        )
        @test_throws DomainError lieb_ando(
            zeros(3, 3),
            zeros(3, 3),
            eye(3),
            5 // 2,
        )
        @test_throws DomainError lieb_ando(
            zeros(3, 3),
            Variable(3, 3),
            eye(3),
            -3 // 2,
        )
        @test_throws DomainError lieb_ando(
            zeros(3, 3),
            Variable(3, 3),
            eye(3),
            5 // 2,
        )
        @test_throws DomainError lieb_ando(
            Variable(3, 3),
            zeros(3, 3),
            eye(3),
            -3 // 2,
        )
        @test_throws DomainError lieb_ando(
            Variable(3, 3),
            zeros(3, 3),
            eye(3),
            5 // 2,
        )
        @test_throws DomainError lieb_ando(
            Variable(3, 3),
            Variable(3, 3),
            eye(3),
            -3 // 2,
        )
        @test_throws DomainError lieb_ando(
            Variable(3, 3),
            Variable(3, 3),
            eye(3),
            5 // 2,
        )
    end

    for n in [2, 3]
        for t in [1 // 2, 1 // 4, 3 // 4, 1 // 8, 3 // 2, 5 // 4]
            for cplx in [false, true]
                # @show n, t, cplx

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
                A = A * A'
                B = B * B'
                A = A / tr(A)
                B = B / tr(B)

                QtAB = lieb_ando(A, B, eye(n), t)

                objective = lieb_ando(X, B, eye(n), t)
                if t >= 0 && t <= 1
                    p = maximize(objective, [X == A]; numeric_type = T)
                else
                    p = minimize(objective, [X == A]; numeric_type = T)
                end
                handle_problem!(p)
                if test
                    if 0 <= t <= 1
                        @test vexity(objective) == ConcaveVexity()
                    else
                        @test vexity(objective) == ConvexVexity()
                    end
                    @test p.optval ≈ QtAB atol = atol * 5 rtol = rtol
                end

                objective = lieb_ando(A, Y, eye(n), t)
                if 0 <= t <= 1
                    p = maximize(objective, [Y == B]; numeric_type = T)
                else
                    p = minimize(objective, [Y == B]; numeric_type = T)
                end
                handle_problem!(p)
                if test
                    if 0 <= t <= 1
                        @test vexity(objective) == ConcaveVexity()
                    else
                        @test vexity(objective) == ConvexVexity()
                    end
                    @test p.optval ≈ QtAB atol = atol * 5 rtol = rtol
                end

                objective = lieb_ando(X, Y, eye(n), t)
                if t >= 0 && t <= 1
                    p = maximize(objective, [X == A, Y == B]; numeric_type = T)
                else
                    p = minimize(objective, [X == A, Y == B]; numeric_type = T)
                end
                handle_problem!(p)
                if test
                    @test p.optval ≈ QtAB atol = atol * 5 rtol = rtol
                end
            end
        end
    end
end

@add_problem sdp function sdp_min_maxeig_canon_lmi(
    handle_problem!,
    ::Val{test},
    atol,
    rtol,
    ::Type{T},
) where {T,test}
    # Minimize the maximum eigenvalue of an affine matrix function of a vector argument. Formulated
    # as a linear optimization subject to a canonical LMI (SDP) constraint:
    #
    # minimize      λ
    # subject to    A(x) ⪯ λ*I
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

    A = A₀ + A₁ * x[1] + A₂ * x[2]

    p = minimize(λ, A ⪯ λ * eye(2); numeric_type = T)

    handle_problem!(p)

    if test
        @test evaluate(λ) ≈ 2.400000051025101 atol = atol rtol = rtol
        @test evaluate(x) ≈ [3.0000000535867315, -0.4000000018585541] atol =
            atol rtol = rtol
        @test evaluate(A) ≈ [
            2.400000046152515 -9.292770553059881e-9
            -9.292770553059881e-9 2.399999957564593
        ] atol = atol rtol = rtol
    end
end

# Using the formulation from:
# https://discourse.julialang.org/t/minimisation-of-operator-norm-solution-does-not-match-evaluated-value-convex-jl/102319/2
@add_problem sdp function sdp_issue_510(
    handle_problem!,
    ::Val{test},
    atol,
    rtol,
    ::Type{T},
) where {T,test}
    d₁ = 2
    d₂ = 2
    H = rand(ComplexF64, d₁ * d₂, d₁ * d₂)
    U = exp(im * π * (H + H'))
    K = rand(ComplexF64, d₁ * d₂, d₁ * d₂)
    K *= d₁ / tr(K' * K)
    ρ = Semidefinite(d₂)
    J = sum(
        kron(
            ((1:d₁) .== j) * ((1:d₁) .== k)',
            partialtrace(
                U' *
                (
                    K *
                    kron(
                        ((1:d₁) .== j) * ((1:d₁) .== k)',
                        LinearAlgebra.I(d₂),
                    ) *
                    K' - kron(((1:d₁) .== j) * ((1:d₁) .== k)', ρ)
                ) *
                U,
                2,
                [d₁, d₂],
            ),
        ) for j in 1:d₁, k in 1:d₁
    )
    constraints = [tr(ρ) == 1]
    p = minimize(LinearAlgebra.opnorm(J, Inf), constraints; numeric_type = T)
    handle_problem!(p)
    if test
        @test p.optval ≈ LinearAlgebra.opnorm(evaluate(J), Inf) atol = atol rtol =
            rtol
        @test tr(evaluate(ρ)) ≈ 1 atol = atol rtol = rtol
    end
end
