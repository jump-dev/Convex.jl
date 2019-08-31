@add_benchmark function affine_negate_atom(handle_problem)
    x = Variable()
    p = minimize(-x, [x <= 0])
    handle_problem(p)
end

@add_benchmark function affine_kron_atom(handle_problem)
    x = ComplexVariable(3, 3)
    y = [1.0 2.0; 3.0 4.0]
    p = satisfy(kron(x,y) == kron(eye(3), y))
    handle_problem(p)
end

@add_benchmark function affine_mult_atom(handle_problem)
    probs=[
        let
            x = Variable(1)
            p = minimize(2.0 * x, [x >= 2, x <= 4])
        end,
        let
            x = Variable(2)
            A = 1.5 * eye(2)
            p = minimize([2 2] * x, [A * x >= [1.1; 1.1]])
        end,
        let
            y = Variable(1)
            x = Variable(3)
            z = [1.0, 2.0, 3.0] * y
            k = -y * [1.0, 2.0, 3.0]
            c = [y <= 3.0, y >= 0.0, x >= ones(3),  k <= x, x <= z]
            o = 3 * y
            p = Problem(:minimize, o, c)
        end,
        let
            x = ComplexVariable(2,2)
            p = minimize( real( [1.0im, 0.0]' * x * [1.0im, 0.0] ), [ x == [1.0 0.0; 0.0 1.0] ])
        end
    ]
    handle_problem.(probs)
end

@add_benchmark function affine_dot_atom(handle_problem)
    probs = [
        let
            x = Variable(2)
            p = minimize(dot([2.0; 2.0], x), x >= [1.1; 1.1])
        end,
        let
            x = Variable(2,2)
            p = minimize(dot(fill(2.0, (2,2)), x), x >= 1.1)
        end
    ]

    handle_problem.(probs)
end


@add_benchmark function affine_add_atom(handle_problem)
    
    probs = [
        let
            x = Variable(1)
            y = Variable(1)
            p = minimize(x + y, [x >= 3, y >= 2])
        end,
        let
            x = Variable(1)
            p = minimize(x, [eye(2) + x >= eye(2)])
        end,
        let
            y = Variable()
            p = minimize(y - 5, y >= -1)
        end ]

    handle_problem.(probs)
end



@add_benchmark function affine_transpose_atom(handle_problem)
    probs = [
        let
            x = Variable(2)
            c = ones(2, 1)
            p = minimize(x' * c, x >= 1)
        end,
        let
            X = Variable(2, 2)
            c = ones(2, 1)
            p = minimize(c' * X' * c, [X >= ones(2, 2)])
        end,
        let
            rows = 2
            cols = 3
            r = rand(rows, cols)
            r_2 = rand(cols, rows)
            x = Variable(rows, cols)
            c = ones(1, cols)
            d = ones(rows, 1)
            p = minimize(c * x' * d + d' * x * c' + (c * x''''' * d)',
                        [x' >= r_2, x >= r, x''' >= r_2, x'' >= r])
        end ]

        handle_problem.(probs)
end

@add_benchmark function affine_index_atom(handle_problem)
        probs = [
        let
            x = Variable(2)
            p = minimize(x[1] + x[2], [x >= 1])
        end,
        let
            x = Variable(3)
            I = [true true false]
            p = minimize(sum(x[I]), [x >= 1])
        end,
        let
            rows = 6
            cols = 8
            n = 2
            X = Variable(rows, cols)
            A = randn(rows, cols)
            c = rand(1, n)
            p = minimize(c * X[1:n, 5:5+n-1]' * c', X >= A)
        end]
        handle_problem.(probs)
end

@add_benchmark function affine_sum_atom(handle_problem)
        probs = [
        let
            x = Variable(2,2)
            p = minimize(sum(x) - 2*x[1,1], x>=1, x[1,1]<=2)
        end,
        let
            x = Variable(2,2)
            p = minimize(sum(x) - 2*x[1,1], x>=1, x[1,1]<=2)
        end
        ]
    return handle_problem.(probs)
end

@add_benchmark function affine_diag_atom(handle_problem)
    probs = [
        let
            x = Variable(2,2)
            p = minimize(sum(diag(x,1)), x >= 1)
        end,
        let
            x = Variable(4, 4)
            p = minimize(sum(diag(x)), x >= 2)
        end,
    ]
    return handle_problem.(probs)

end

@add_benchmark function affine_tr_atom(handle_problem)
    x = Variable(2,2)
    p = minimize(tr(x), x >= 1)
    handle_problem(p)
end
