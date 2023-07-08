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
            add_constraints_to_context(a >= 0, context)
        elseif sign(a) == Negative()
            add_constraints_to_context(a <= 0, context)
        end

        if vartype(a) == BinVar
            MOI.add_constraints(
                context.model,
                [MOI.SingleVariable(i) for i in var_inds],
                [MOI.ZeroOne() for i in var_inds],
            )
        elseif vartype(a) == IntVar
            MOI.add_constraints(
                context.model,
                [MOI.SingleVariable(i) for i in var_inds],
                [MOI.Integer() for i in var_inds],
            )
        end

        for constraint in constraints(a)
            add_constraints_to_context(constraint, context)
        end
    end
    return MOI.VectorOfVariables(var_inds)
end

function to_tape(v::MOI.VectorOfVariables, context::Context{T}) where {T}
    var_inds = v.variables
    d = length(var_inds)
    if USE_SPARSE2()
        return SparseVAFTape(
            [SparseAffineOperation(sparse(one(T) * I, d, d), zeros(T, d))],
            var_inds,
        )
    else
        return VAFTape(
            tuple(AffineOperation(one(T) * I, zeros(T, d))),
            var_inds,
        )
    end
end

# get the usual tape
function template(a::AbstractVariable, context::Context)
    if vexity(a) == ConstVexity()
        return template(constant(evaluate(a)), context)
    end
    return to_tape(_template(a, context), context)
end

function template(c::ComplexVariable, context::Context)
    if vexity(c) == ConstVexity()
        return template(constant(evaluate(c)), context)
    end
    re = template(c.real_var, context)
    im = template(c.imag_var, context)
    for constraint in constraints(c)
        add_constraints_to_context(constraint, context)
    end
    return ComplexTape(re, im)
end