### (Primal) exponential cone constraint ExpConstraint(x,y,z) => y exp(x/y) <= z & y>=0
mutable struct ExpConstraint <: Constraint
    head::Symbol
    id_hash::UInt64
    children::Tuple{AbstractExpr,AbstractExpr,AbstractExpr} # (x, y, z)
    size::Tuple{Int,Int}
    dual::ValueOrNothing

    function ExpConstraint(x::AbstractExpr, y::AbstractExpr, z::AbstractExpr)
        @assert(
            x.size == y.size == z.size,
            "Exponential constraint requires x, y, and z to be of same size"
        )
        # @assert(x.size == (1,1),
        #         "Exponential constraint requires x, y, and z to be scalar for now")
        sz = x.size
        id_hash = hash((x, y, z, :exp))
        return new(:exp, id_hash, (x, y, z), sz, nothing)
    end
end

ExpConstraint(x::AbstractExpr, y, z::AbstractExpr) = ExpConstraint(x, constant(y), z)
ExpConstraint(x::AbstractExpr, y::AbstractExpr, z) = ExpConstraint(x, y, constant(z))
ExpConstraint(x, y::AbstractExpr, z::AbstractExpr) = ExpConstraint(constant(x), y, z)

function vexity(c::ExpConstraint)
    # TODO: check these...
    if vexity(c.x) == ConcaveVexity()
        error("Exponential constraint requires x to be convex")
    end
    if vexity(c.y) != ConstVexity()
        error("Exponential constraint requires y to be constant")
    end
    if vexity(c.z) == ConvexVexity()
        error("Exponential constraint requires z to be concave")
    end
    return ConvexVexity()
end

function _add_constraints_to_context(c::ExpConstraint, context::Context{T}) where {T}
    x, y, z = c.children
    t = a -> template(a, context)
    inds = [
        begin
            terms = (t(x[i, j]), t(y[i, j]), t(z[i, j]))
            obj = operate(vcat, T, terms...)
            MOI_add_constraint(context.model, obj, MOI.ExponentialCone())
        end for i = 1:size(x, 1), j = 1:size(x, 2)
    ]
    context.constr_to_moi_inds[c] = inds
    return nothing
end

function populate_dual!(model::MOI.ModelLike, constr::ExpConstraint, MOI_constr_indices)
    x = first(children(c))
    constr.dual = output([ MOI.get(model, MOI.ConstraintDual(), MOI_constr_indices[i,j]) for i = 1:size(x, 1), j = 1:size(x, 2)])
    return nothing
end
