@add_benchmark function sdp_nuclear_norm_atom(handle_problem)
    y = Semidefinite(3)
    p = minimize(nuclearnorm(y), y[2,1]<=4, y[2,2]>=3, y[3,3]<=2)
    handle_problem(p)
end

@add_benchmark function sdp_operator_norm_atom(handle_problem)
    y = Variable((3,3))
    p = minimize(opnorm(y), y[2,1]<=4, y[2,2]>=3, sum(y)>=12)
    handle_problem(p)
end

@add_benchmark function sdp_sigma_max_atom(handle_problem)
    y = Variable((3,3))
    p = minimize(sigmamax(y), y[2,1]<=4, y[2,2]>=3, sum(y)>=12)
    handle_problem(p)
end

@add_benchmark function sdp_lambda_max_atom(handle_problem)
    y = Semidefinite(3)
    p = minimize(lambdamax(y), y[1,1]>=4)
    handle_problem(p)
end

@add_benchmark function sdp_lambda_min_atom(handle_problem)
    y = Semidefinite(3)
    p = maximize(lambdamin(y), tr(y)<=6)
    handle_problem(p)
end


@add_benchmark function sdp_matrix_frac_atom(handle_problem)
    probs = [
        let
            x = [1, 2, 3]
            P = Variable(3, 3)
            p = minimize(matrixfrac(x, P), P <= 2*eye(3), P >= 0.5 * eye(3))
        end,
        let
            x = Variable(3)
            P = Variable(3, 3)
            p = minimize(matrixfrac(x, P), lambdamax(P) <= 2, x[1] >= 1)
        end
    ]

    handle_problem.(probs)
end

@add_benchmark function sdp_sum_squares_atom(handle_problem)
    n = 10
    A = rand(n,n) + im*rand(n,n)
    A = A + A' # now A is hermitian
    x = ComplexVariable(n,n)
    objective = sumsquares(A - x)
    c1 = x in :SDP
    p = minimize(objective, c1)
    handle_problem(p)
end

