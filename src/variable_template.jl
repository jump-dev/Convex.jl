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

    return to_vov(var_inds)
end

# Real case
to_vov(var_inds::Vector{MOI.VariableIndex}) = MOI.VectorOfVariables(var_inds)

# Complex case
function to_vov(
    (v1, v2)::Tuple{Vector{MOI.VariableIndex},Vector{MOI.VariableIndex}},
)
    return (to_vov(v1), to_vov(v2))
end

function to_tape(v::MOI.VectorOfVariables, ::Context{T}) where {T}
    var_inds = v.variables
    d = length(var_inds)
    return SparseTape(
        [SparseAffineOperation(gbidentity(T, d), GBVector{T,T}(d))],
        var_inds,
    )
end

# Complex case
function to_tape(
    (v1, v2)::Tuple{MOI.VectorOfVariables,MOI.VectorOfVariables},
    context::Context{T},
) where {T}
    return ComplexTape(to_tape(v1, context), to_tape(v2, context))
end

# get the usual tape
function _conic_form!(context::Context, a::AbstractVariable)
    if vexity(a) == ConstVexity()
        return conic_form!(context, constant(evaluate(a)))
    end
    return to_tape(_template(a, context), context)
end

function conic_form!(context::Context, a::AbstractExpr)

    # Nicer implementation
    # d = context.conic_form_cache
    # return get!(() -> _conic_form!(context, a), d, a)

    # Avoid closure
    wkh = context.conic_form_cache
    default = () -> _conic_form!(context, a)
    key = a
    return Base.@lock wkh.lock begin
        get!(default, wkh.ht, DataStructures.WeakRefForWeakDict(key))
    end

    return
end
