USE_SPARSE2() = true

# It might be useful to get a direct VOV sometimes...
function _template(a::AbstractVariable, context::Context{T}) where {T}
    first_cache = false
    var_inds = get!(context.var_id_to_moi_indices, a.id_hash) do
        first_cache = true
        return add_variables!(context.model, a)
    end

    context.id_to_variables[a.id_hash] = a

    # we only want this to run once, when the variable is first added,
    # and after `var_id_to_moi_indices` is populated
    if first_cache
        if sign(a) == Positive()
            add_constraint!(context, a >= 0)
        elseif sign(a) == Negative()
            add_constraint!(context, a <= 0)
        end

        if vartype(a) == BinVar
            MOI.add_constraints(
                context.model,
                var_inds,
                [MOI.ZeroOne() for i in var_inds],
            )
        elseif vartype(a) == IntVar
            MOI.add_constraints(
                context.model,
                var_inds,
                [MOI.Integer() for i in var_inds],
            )
        end

        for constraint in constraints(a)
            add_constraint!(context, constraint)
        end
    end
    return MOI.VectorOfVariables(var_inds)
end

function to_tape(v::MOI.VectorOfVariables, context::Context{T}) where {T}
    var_inds = v.variables
    d = length(var_inds)
    return SparseVAFTape(
        [SparseAffineOperation(sparse(one(T) * I, d, d), zeros(T, d))],
        var_inds,
    )
end

# get the usual tape
function conic_form!(context::Context, a::AbstractVariable)
    if vexity(a) == ConstVexity()
        return conic_form!(context, constant(evaluate(a)))
    end
    return to_tape(_template(a, context), context)
end

function conic_form!(context::Context, c::ComplexVariable)
    if vexity(c) == ConstVexity()
        return conic_form!(context, constant(evaluate(c)))
    end
    re = conic_form!(context, c.real_var)
    im = conic_form!(context, c.imag_var)
    for constraint in constraints(c)
        add_constraint!(context, constraint)
    end
    return ComplexTape(re, im)
end
