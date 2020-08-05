@add_problem affine function affine_negate_atom(handle_problem!, ::Val{test}, atol, rtol, ::Type{T}) where {T, test}
    x = Variable()
    p = minimize(-x, [x <= 0]; numeric_type = T)

    if test
        @test vexity(p) == AffineVexity()
    end
    handle_problem!(p)
    if test
        @test p.optval ≈ 0 atol=atol rtol=rtol
        @test evaluate(-x) ≈ 0 atol=atol rtol=rtol
    end
end

@add_problem affine function affine_kron_atom(handle_problem!, ::Val{test}, atol, rtol, ::Type{T}) where {T, test}
    x = ComplexVariable(3, 3)
    y = [1.0 2.0; 3.0 4.0]
    if test
        @test size(kron(x, y)) == (6, 6)
        @test size(kron(y, x)) == (6, 6)
    end
end

@add_problem affine function affine_multiply_atom(handle_problem!, ::Val{test}, atol, rtol, ::Type{T}) where {T, test}
    x = Variable(1)
    p = minimize(2.0 * x, [x >= 2, x <= 4]; numeric_type = T)

    if test
        @test vexity(p) == AffineVexity()
    end
    handle_problem!(p)
    if test
        @test p.optval ≈ 4 atol=atol rtol=rtol
        @test (evaluate(2.0x))[1] ≈ 4 atol=atol rtol=rtol
    end

    x = Variable(2)
    A = 1.5 * eye(2)
    p = minimize([2 2] * x, [A * x >= [1.1; 1.1]]; numeric_type = T)

    if test
        @test vexity(p) == AffineVexity()
    end
    handle_problem!(p)
    if test
        @test p.optval ≈ 2.93333 atol=atol rtol=rtol
        @test (evaluate([2 2] * x))[1] ≈ 2.93333 atol=atol rtol=rtol
        @test vec(evaluate(A * x)) ≈ [1.1; 1.1] atol=atol rtol=rtol
    end

    y = Variable(1)
    x = Variable(3)
    z = [1.0, 2.0, 3.0] * y
    k = -y * [1.0, 2.0, 3.0]
    c = [y <= 3.0, y >= 0.0, x >= ones(3), k <= x, x <= z]
    o = 3 * y
    p = Problem{T}(:minimize, o, c)
    if test
        @test vexity(p) == AffineVexity()
    end
    handle_problem!(p)
    if test
        @test p.optval ≈ 3 atol=atol rtol=rtol
    end

    p = Problem{T}(:minimize, o, c...)
    if test
        @test vexity(p) == AffineVexity()
    end
    handle_problem!(p)
    if test
        @test p.optval ≈ 3 atol=atol rtol=rtol
    end

    # Check #274
    x = ComplexVariable(2,2)
    p = minimize( real( [1.0im, 0.0]' * x * [1.0im, 0.0] ), [ x == [1.0 0.0; 0.0 1.0] ]; numeric_type = T)

    handle_problem!(p)
    if test
        @test p.optval ≈ 1.0 atol=atol rtol=rtol
    end
end

@add_problem affine function affine_dot_atom(handle_problem!, ::Val{test}, atol, rtol, ::Type{T}) where {T, test}
    x = Variable(2)
    p = minimize(dot([2.0; 2.0], x), x >= [1.1; 1.1]; numeric_type = T)

    if test
        @test vexity(p) == AffineVexity()
    end
    handle_problem!(p)
    if test
        @test p.optval ≈ 4.4 atol=atol rtol=rtol
        @test (evaluate(dot([2.0; 2.0], x)))[1] ≈ 4.4 atol=atol rtol=rtol
    end
end

@add_problem affine function affine_dot_atom_for_matrix_variables(handle_problem!, ::Val{test}, atol, rtol, ::Type{T}) where {T, test}
    x = Variable(2,2)
    p = minimize(dot(fill(2.0, (2,2)), x), x >= 1.1; numeric_type = T)

    if test
        @test vexity(p) == AffineVexity()
    end
    handle_problem!(p)
    if test
        @test p.optval ≈ 8.8 atol=atol rtol=rtol
        @test (evaluate(dot(fill(2.0, (2, 2)), x)))[1] ≈ 8.8 atol=atol rtol=rtol
    end
end

@add_problem affine function affine_add_atom(handle_problem!, ::Val{test}, atol, rtol, ::Type{T}) where {T, test}
    x = Variable(1)
    y = Variable(1)
    p = minimize(x + y, [x >= 3, y >= 2]; numeric_type = T)

    if test
        @test vexity(p) == AffineVexity()
    end
    handle_problem!(p)
    if test
        @test p.optval ≈ 5 atol=atol rtol=rtol
        @test evaluate(x + y) ≈ 5 atol=atol rtol=rtol
    end

    x = Variable(1)
    p = minimize(x, [eye(2) + x*ones(2,2) >= eye(2)]; numeric_type = T)

    if test
        @test vexity(p) == AffineVexity()
    end
    handle_problem!(p)
    if test
        @test p.optval ≈ 0 atol=atol rtol=rtol
        @test evaluate(eye(2) + x*ones(2,2)) ≈ eye(2) atol=atol rtol=rtol
    end

    y = Variable()
    p = minimize(y - 5, y >= -1; numeric_type = T)

    if test
        @test vexity(p) == AffineVexity()
    end
    handle_problem!(p)
    if test
        @test p.optval ≈ -6 atol=atol rtol=rtol
        @test evaluate(y - 5) ≈ -6 atol=atol rtol=rtol
    end
end

@add_problem affine function affine_transpose_atom(handle_problem!, ::Val{test}, atol, rtol, ::Type{T}) where {T, test}
    x = Variable(2)
    c = ones(2, 1)
    p = minimize(x' * c, x >= 1; numeric_type = T)

    if test
        @test vexity(p) == AffineVexity()
    end
    handle_problem!(p)
    if test
        @test p.optval ≈ 2 atol=atol rtol=rtol
        @test (evaluate(x' * c))[1] ≈ 2 atol=atol rtol=rtol
    end

    X = Variable(2, 2)
    c = ones(2, 1)
    p = minimize(c' * X' * c, [X >= ones(2, 2)]; numeric_type = T)

    if test
        @test vexity(p) == AffineVexity()
    end
    handle_problem!(p)
    if test
        @test p.optval ≈ 4 atol=atol rtol=rtol
        @test (evaluate(c' * X' * c))[1] ≈ 4 atol=atol rtol=rtol
    end

    rows = 2
    cols = 3
    r = rand(rows, cols)
    r_2 = rand(cols, rows)
    x = Variable(rows, cols)
    c = ones(1, cols)
    d = ones(rows, 1)
    p = minimize(c * x' * d + d' * x * c' + (c * x''''' * d)',
                [x' >= r_2, x >= r, x''' >= r_2, x'' >= r]; numeric_type = T)
    if test
        @test vexity(p) == AffineVexity()
    end
    handle_problem!(p)
    s = sum(max.(r, r_2')) * 3
    if test
        @test p.optval ≈ s atol=atol rtol=rtol
        @test (evaluate(c * x' * d + d' * x * c' + (c * ((((x')')')')' * d)'))[1] ≈ s atol=atol rtol=rtol
    end
end

@add_problem affine function affine_index_atom(handle_problem!, ::Val{test}, atol, rtol, ::Type{T}) where {T, test}
    x = Variable(2)
    p = minimize(x[1] + x[2], [x >= 1]; numeric_type = T)

    if test
        @test vexity(p) == AffineVexity()
    end
    handle_problem!(p)
    if test
        @test p.optval ≈ 2 atol=atol rtol=rtol
        @test (evaluate(x[1] + x[2]))[1] ≈ 2 atol=atol rtol=rtol
    end

    x = Variable(3)
    I = [true true false]
    p = minimize(sum(x[I]), [x >= 1]; numeric_type = T)

    if test
        @test vexity(p) == AffineVexity()
    end
    handle_problem!(p)
    if test
        @test p.optval ≈ 2 atol=atol rtol=rtol
        @test (evaluate(sum(x[I])))[1] ≈ 2 atol=atol rtol=rtol
    end

    rows = 6
    cols = 8
    n = 2
    X = Variable(rows, cols)
    A = randn(rows, cols)
    c = rand(1, n)
    p = minimize(c * X[1:n, 5:5+n-1]' * c', X >= A; numeric_type = T)

    if test
        @test vexity(p) == AffineVexity()
    end
    handle_problem!(p)
    s = c * A[1:n, 5:5+n-1]' * c'
    if test
        @test p.optval ≈ s[1] atol=atol rtol=rtol
        @test evaluate(c * (X[1:n, 5:(5 + n) - 1])' * c') ≈ s atol=atol rtol=rtol
    end
end

@add_problem affine function affine_sum_atom(handle_problem!, ::Val{test}, atol, rtol, ::Type{T}) where {T, test}
    x = Variable(2,2)
    p = minimize(sum(x), x>=1; numeric_type = T)

    if test
        @test vexity(p) == AffineVexity()
    end
    handle_problem!(p)
    if test
        @test p.optval ≈ 4 atol=atol rtol=rtol
        @test evaluate(sum(x)) ≈ 4 atol=atol rtol=rtol
    end

    x = Variable(2,2)
    p = minimize(sum(x) - 2*x[1,1], x>=1, x[1,1]<=2; numeric_type = T)

    if test
        @test vexity(p) == AffineVexity()
    end
    handle_problem!(p)
    if test
        @test p.optval ≈ 1 atol=atol rtol=rtol
        @test (evaluate(sum(x) - 2 * x[1, 1]))[1] ≈ 1 atol=atol rtol=rtol
    end

    x = Variable(10)
    a = rand(10, 1)
    p = maximize(sum(x[2:6]), x <= a; numeric_type = T)

    if test
        @test vexity(p) == AffineVexity()
    end
    handle_problem!(p)
    if test
        @test p.optval ≈ sum(a[2:6]) atol=atol rtol=rtol
        @test evaluate(sum(x[2:6])) ≈ sum(a[2:6]) atol=atol rtol=rtol
    end
end

@add_problem affine function affine_diag_atom(handle_problem!, ::Val{test}, atol, rtol, ::Type{T}) where {T, test}
    x = Variable(2,2)
    p = minimize(sum(diag(x,1)), x >= 1; numeric_type = T)

    if test
        @test vexity(p) == AffineVexity()
    end
    handle_problem!(p)
    if test
        @test p.optval ≈ 1 atol=atol rtol=rtol
        @test evaluate(sum(diag(x, 1))) ≈ 1 atol=atol rtol=rtol
    end

    x = Variable(4, 4)
    p = minimize(sum(diag(x)), x >= 2; numeric_type = T)

    if test
        @test vexity(p) == AffineVexity()
    end
    handle_problem!(p)
    if test
        @test p.optval ≈ 8 atol=atol rtol=rtol
        @test evaluate(sum(diag(x))) ≈ 8 atol=atol rtol=rtol
    end
end

@add_problem affine function affine_trace_atom(handle_problem!, ::Val{test}, atol, rtol, ::Type{T}) where {T, test}
    x = Variable(2,2)
    p = minimize(tr(x), x >= 1; numeric_type = T)

    if test
        @test vexity(p) == AffineVexity()
    end
    handle_problem!(p)
    if test
        @test p.optval ≈ 2 atol=atol rtol=rtol
        @test evaluate(tr(x)) ≈ 2 atol=atol rtol=rtol
    end
end

@add_problem affine function affine_dot_multiply_atom(handle_problem!, ::Val{test}, atol, rtol, ::Type{T}) where {T, test}
    x = Variable(3)
    p = maximize(sum(dot(*)(x,[1,2,3])), x<=1; numeric_type = T)

    if test
        @test vexity(p) == AffineVexity()
    end
    handle_problem!(p)
    if test
        @test p.optval ≈ 6 atol=atol rtol=rtol
        @test evaluate(sum((dot(*))(x, [1, 2, 3]))) ≈ 6 atol=atol rtol=rtol
    end

    x = Variable(3, 3)
    p = maximize(sum(dot(*)(x,eye(3))), x<=1; numeric_type = T)

    if test
        @test vexity(p) == AffineVexity()
    end
    handle_problem!(p)
    if test
        @test p.optval ≈ 3 atol=atol rtol=rtol
        @test evaluate(sum((dot(*))(x, eye(3)))) ≈ 3 atol=atol rtol=rtol
    end

    x = Variable(5, 5)
    p = minimize(x[1, 1], dot(*)(3,x) >= 3; numeric_type = T)

    if test
        @test vexity(p) == AffineVexity()
    end
    handle_problem!(p)
    if test
        @test p.optval ≈ 1 atol=atol rtol=rtol
        @test (evaluate(x[1, 1]))[1] ≈ 1 atol=atol rtol=rtol
    end

    x = Variable(3,1)
    p = minimize(sum(dot(*)(ones(3,3), x)), x>=1; numeric_type = T)

    if test
        @test vexity(p) == AffineVexity()
    end
    handle_problem!(p)
    if test
        @test p.optval ≈ 9 atol=atol rtol=rtol
        @test (evaluate(x[1, 1]))[1] ≈ 1 atol=atol rtol=rtol
    end

    x = Variable(1,3)
    p = minimize(sum(dot(*)(ones(3,3), x)), x>=1; numeric_type = T)

    if test
        @test vexity(p) == AffineVexity()
    end
    handle_problem!(p)
    if test
        @test p.optval ≈ 9 atol=atol rtol=rtol
        @test (evaluate(x[1, 1]))[1] ≈ 1 atol=atol rtol=rtol
    end

    x = Variable(1, 3, Positive())
    p = maximize(sum(dot(/)(x,[1 2 3])), x<=1; numeric_type = T)

    if test
        @test vexity(p) == AffineVexity()
    end
    handle_problem!(p)
    if test
        @test p.optval ≈ 11 / 6 atol=atol rtol=rtol
        @test evaluate(sum((dot(/))(x, [1 2 3]))) ≈ 11 / 6 atol=atol rtol=rtol
    end

    # Broadcast fusion works
    x = Variable(5, 5)
    a = 2.0 .* x .* ones(Int, 5)

    xval = rand(5,5)
    if test
        @test a isa Convex.DotMultiplyAtom
    end
end

@add_problem affine function affine_reshape_atom(handle_problem!, ::Val{test}, atol, rtol, ::Type{T}) where {T, test}
    A = rand(2, 3)
    X = Variable(3, 2)
    c = rand()
    p = minimize(sum(reshape(X, 2, 3) + A), X >= c; numeric_type = T)

    if test
        @test vexity(p) == AffineVexity()
    end
    handle_problem!(p)
    if test
        @test p.optval ≈ sum(A .+ c) atol=atol rtol=rtol
        @test evaluate(sum(reshape(X, 2, 3) + A)) ≈ sum(A .+ c) atol=atol rtol=rtol
    end

    b = rand(6)
    p = minimize(sum(vec(X) + b), X >= c; numeric_type = T)

    if test
        @test vexity(p) == AffineVexity()
    end
    handle_problem!(p)
    if test
        @test p.optval ≈ sum(b .+ c) atol=atol rtol=rtol
        @test evaluate(sum(vec(X) + b)) ≈ sum(b .+ c) atol=atol rtol=rtol
    end

    x = Variable(4, 4)
    c = ones(16, 1)
    reshaped = reshape(x, 16, 1)
    a = collect(1:16)
    p = minimize(c' * reshaped, reshaped >= a; numeric_type = T)

    if test
        @test vexity(p) == AffineVexity()
    end
    handle_problem!(p)
    # TODO: why is accuracy lower here?
    if test
        @test p.optval ≈ 136 atol=10atol atol=atol rtol=rtol
        @test (evaluate(c' * reshaped))[1] ≈ 136 atol=10atol atol=atol rtol=rtol
    end
end

@add_problem affine function affine_hcat_atom(handle_problem!, ::Val{test}, atol, rtol, ::Type{T}) where {T, test}
    x = Variable(4, 4)
    y = Variable(4, 6)
    p = maximize(sum(x) + sum([y fill(4.0, 4)]), [x y fill(2.0, (4, 2))] <= 2; numeric_type = T)

    if test
        @test vexity(p) == AffineVexity()
    end
    handle_problem!(p)
    if test
        @test p.optval ≈ 96 atol=atol rtol=rtol
        @test evaluate(sum(x) + sum([y fill(4.0, 4)])) ≈ 96 atol=atol rtol=rtol
        @test evaluate([x y fill(2.0, (4, 2))]) ≈ fill(2.0, (4, 12)) atol=atol rtol=rtol
    end
end

@add_problem affine function affine_vcat_atom(handle_problem!, ::Val{test}, atol, rtol, ::Type{T}) where {T, test}
    x = Variable(4, 4)
    y = Variable(4, 6)

    # TODO: fix dimension mismatch [y 4*eye(4); x -ones(4, 6)]
    p = maximize(sum(x) + sum([y 4*eye(4); x -ones(4, 6)]), [x;y'] <= 2; numeric_type = T)

    if test
        @test vexity(p) == AffineVexity()
    end
    handle_problem!(p)
    # TODO: why is accuracy lower here?
    if test
        @test p.optval ≈ 104 atol=10atol atol=atol rtol=rtol
        @test evaluate(sum(x) + sum([y 4 * eye(4); x -(ones(4, 6))])) ≈ 104 atol=10atol atol=atol rtol=rtol
        @test evaluate([x; y']) ≈ 2 * ones(10, 4) atol=atol rtol=rtol
    end
end

@add_problem affine function affine_Diagonal_atom(handle_problem!, ::Val{test}, atol, rtol, ::Type{T}) where {T, test}
    x = Variable(2, 2)
    if test
        @test_throws ArgumentError Diagonal(x)
    end

    x = Variable(4)
    p = minimize(sum(Diagonal(x)), x == [1, 2, 3, 4]; numeric_type = T)

    if test
        @test vexity(p) == AffineVexity()
    end
    handle_problem!(p)
    if test
        @test p.optval ≈ 10 atol=atol rtol=rtol
        @test all(abs.(evaluate(Diagonal(x)) - Diagonal([1, 2, 3, 4])) .<= atol)
    end

    x = Variable(3)
    c = [1; 2; 3]
    p = minimize(c' * Diagonal(x) * c, x >= 1, sum(x) == 10; numeric_type = T)

    if test
        @test vexity(p) == AffineVexity()
    end
    handle_problem!(p)
    if test
        @test p.optval ≈ 21 atol=atol rtol=rtol
    end

    x = Variable(3)
    p = minimize(sum(x), x >= 1, Diagonal(x)[1, 2] == 1; numeric_type = T)

    handle_problem!(p)
    if test
        @test p.status != MOI.OPTIMAL
    end
end

@add_problem affine function affine_conv_atom(handle_problem!, ::Val{test}, atol, rtol, ::Type{T}) where {T, test}
    x = Variable(3)
    h = [1, -1]
    p = minimize(sum(conv(h, x)) + sum(x), x >= 1, x <= 2; numeric_type = T)

    if test
        @test vexity(p) == AffineVexity()
    end
    handle_problem!(p)
    if test
        @test p.optval ≈ 3 atol=atol rtol=rtol
        @test evaluate(sum(conv(h, x))) ≈ 0 atol=atol rtol=rtol
    end

    x = Variable(3)
    h = [1, -1]
    p = minimize(sum(conv(x, h)) + sum(x), x >= 1, x <= 2; numeric_type = T)

    if test
        @test vexity(p) == AffineVexity()
    end
    handle_problem!(p)
    if test
        @test p.optval ≈ 3 atol=atol rtol=rtol
        @test evaluate(sum(conv(h, x))) ≈ 0 atol=atol rtol=rtol
    end

end

@add_problem affine function affine_satisfy_problems(handle_problem!, ::Val{test}, atol, rtol, ::Type{T}) where {T, test}
    x = Variable()
    p = satisfy(x >= 0; numeric_type = T)

    add_constraints!(p, x >= 1)
    add_constraints!(p, [x >= -1, x <= 4])
    handle_problem!(p)
    if test
        @test p.status == MOI.OPTIMAL
    end

    p = satisfy([x >= 0, x >= 1, x <= 2]; numeric_type = T)

    handle_problem!(p)
    if test
        @test p.status == MOI.OPTIMAL
    end

    p = satisfy([x >= 1, x <= 2]; numeric_type = T)

    handle_problem!(p)
    if test
        @test p.status == MOI.OPTIMAL
    end

    constr = x >= 0
    constr += x >= 1
    constr += x <= 10
    constr2 = x >= 0
    constr2 += [x >= 2, x <= 3] + constr
    p = satisfy(constr; numeric_type = T)

    handle_problem!(p)
    if test
        @test p.status == MOI.OPTIMAL
    end
end

@add_problem affine function affine_dualvalue(handle_problem!, ::Val{test}, atol, rtol, ::Type{T}) where {T, test}
    x = Variable()
    p = minimize(x, x >= 0; numeric_type = T)

    handle_problem!(p)
    if test
        @test p.constraints[1].dual ≈ 1 atol=atol rtol=rtol
    end

    x = Variable()
    p = maximize(x, x <= 0; numeric_type = T)

    handle_problem!(p)
    if test
        @test p.constraints[1].dual ≈ 1 atol=atol rtol=rtol
    end

    x = Variable()
    p = minimize(x, x >= 0, x == 2; numeric_type = T)

    handle_problem!(p)
    if test
        @test p.constraints[1].dual ≈ 0 atol=atol rtol=rtol
        @test abs.(p.constraints[2].dual) ≈ 1 atol=atol rtol=rtol
    end

    x = Variable(2)
    A = 1.5 * eye(2)
    p = minimize(dot([2.0; 2.0], x), [A * x >= [1.1; 1.1]]; numeric_type = T)

    handle_problem!(p)
    if test
        dual = [4/3; 4/3]
        @test all(abs.(p.constraints[1].dual - dual) .<= atol)
    end
end

@add_problem affine function affine_Partial_transpose(handle_problem!, ::Val{test}, atol, rtol, ::Type{T}) where {T, test}
    dims = [2,3,4]
    d = prod(dims)
    A = rand(ComplexF64,2,2)
    B = rand(ComplexF64,3,3)
    C = rand(ComplexF64,4,4)
    M = kron(A,B,C)
    Mt1 = kron(transpose(A),B,C)
    Mt2 = kron(A,transpose(B),C)
    Mt3 = kron(A,B,transpose(C))

    Rt1 = ComplexVariable(d,d)
    Rt2 = ComplexVariable(d,d)
    Rt3 = ComplexVariable(d,d)
    S = rand(ComplexF64,d,d)
    handle_problem!(satisfy(partialtranspose(Rt1, 1, dims) == S; numeric_type = T))

    handle_problem!(satisfy(partialtranspose(Rt2, 2, dims) == S; numeric_type = T))

    handle_problem!(satisfy(partialtranspose(Rt3, 3, dims) == S; numeric_type = T))


        
    if test
        @test partialtranspose(M,1,dims) ≈ Mt1 atol=atol rtol=rtol
        @test partialtranspose(M,2,dims) ≈ Mt2 atol=atol rtol=rtol
        @test partialtranspose(M,3,dims) ≈ Mt3 atol=atol rtol=rtol
        @test partialtranspose(S,1,dims) ≈ evaluate(Rt1) atol=atol rtol=rtol
        @test partialtranspose(S,2,dims) ≈ evaluate(Rt2) atol=atol rtol=rtol
        @test partialtranspose(S,3,dims) ≈ evaluate(Rt3) atol=atol rtol=rtol
    end

    if test
        @test_throws ArgumentError partialtrace(rand(6, 6), 3, [2, 3])
        @test_throws ArgumentError partialtrace(rand(6, 6), 1, [2, 4])
        @test_throws ArgumentError partialtrace(rand(3, 4), 1, [2, 3])
    end
end

@add_problem affine function affine_permuteddims_matrix(handle_problem!, ::Val{test}, atol, rtol, ::Type{T}) where {T, test}
#this function is used in the partial transpose 
    for n in (2, 3, 4, 5)
        dims = ntuple( i -> rand(2:5), n)
        d = prod(dims)
        v = rand(d)
        p = randperm(n)
        out1 = vec(permutedims(reshape(v, dims), p))
        out2 = Convex.permutedims_matrix(dims, p) * v
        if test
            @test out1 ≈ out2 atol=atol rtol=rtol
        end
    end
end
