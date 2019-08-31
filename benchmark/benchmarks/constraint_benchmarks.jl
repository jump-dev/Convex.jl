@add_benchmark function LT_constraints(handle_problem)
    p = satisfy()
    x = [Variable() for _ = 1:1000]
    for (i, xi) in enumerate(x)
        push!(p.constraints, xi <= 1.0 * i)
    end
    return handle_problem(p)
end

@add_benchmark function LT_constraint(handle_problem)
    p = satisfy()
    x = Variable(1000)
    push!(p.constraints, x <= collect(1.0:1000.0))
    return handle_problem(p)
end

@add_benchmark function GT_constraints(handle_problem)
    p = satisfy()
    x = [Variable() for _ = 1:1000]
    for (i, xi) in enumerate(x)
        push!(p.constraints, xi >= 1.0 * i)
    end
    return handle_problem(p)
end

@add_benchmark function GT_constraint(handle_problem)
    p = satisfy()
    x = Variable(1000)
    push!(p.constraints, x >= collect(1.0:1000.0))
    return handle_problem(p)
end


@add_benchmark function equality_constraints(handle_problem)
    p = satisfy()
    x = [Variable() for _ = 1:1000]
    for (i, xi) in enumerate(x)
        push!(p.constraints, xi == 1.0 * i)
    end
    return handle_problem(p)
end

@add_benchmark function equality_constraint(handle_problem)
    p = satisfy()
    x = Variable(1000)
    push!(p.constraints, x == collect(1.0:1000.0))
    return handle_problem(p)
end

@add_benchmark function SDP_constraint(handle_problem)
    p = satisfy()
    x = Variable(44, 44) # 990 vectorized entries
    push!(p.constraints, x ⪰ 0)
    return handle_problem(p)
end

@add_benchmark function SDP_constraints(handle_problem)
    p = satisfy()
    x = [Variable(4, 4) for _ = 1:100] # 1000 total vectorized entries
    for v in x
        push!(p.constraints, v ⪰ 0)
    end
    return handle_problem(p)
end
