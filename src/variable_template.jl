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

        for constraint in get_constraints(a)
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
        SparseAffineOperation(spidentity(T, d), spzeros(T, d)),
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

"""
    new_conic_form!(context::Context, a::AbstractExpr)

Create a new conic form for `a` and return it, assuming that no conic form
for `a` has already been created, that is `!haskey(context, a)` as this is
already checked in [`conic_form!`](@ref) which calls this function.
"""
function new_conic_form! end

# get the usual tape
function new_conic_form!(context::Context, a::AbstractVariable)
    if vexity(a) == ConstVexity()
        return conic_form!(context, constant(evaluate(a)))
    end
    return to_tape(_template(a, context), context)
end

"""
    conic_form!(context::Context, a::AbstractExpr)

Return the conic form for `a`. If it as already been created, it is directly
accessed in `context[a]`, otherwise, it is created by calling
[`Convex.new_conic_form!`](@ref) and then cached in `context` so that the next call
with the same expression does not create a duplicate one.
"""
function conic_form!(context::Context, a::AbstractExpr)

    # Nicer implementation
    d = context.conic_form_cache
    return get!(() -> new_conic_form!(context, a), d, a)

    # Avoid closure
    # wkh = context.conic_form_cache
    # default = () -> new_conic_form!(context, a)
    # key = a
    # return Base.@lock wkh.lock begin
    #     get!(default, wkh.ht, DataStructures.WeakRefForWeakDict(key))
    # end

    return
end
