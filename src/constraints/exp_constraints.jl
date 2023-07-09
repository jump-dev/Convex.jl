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

function ExpConstraint(x::AbstractExpr, y, z::AbstractExpr)
    return ExpConstraint(x, constant(y), z)
end
function ExpConstraint(x::AbstractExpr, y::AbstractExpr, z)
    return ExpConstraint(x, y, constant(z))
end
function ExpConstraint(x, y::AbstractExpr, z::AbstractExpr)
    return ExpConstraint(constant(x), y, z)
end

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

function _add_constraint!(context::Context{T}, c::ExpConstraint) where {T}
    x, y, z = c.children
    t = a -> conic_form!(context, a)
    inds = [
        begin
            terms = (t(x[i, j]), t(y[i, j]), t(z[i, j]))
            obj = operate(vcat, T, terms...)
            MOI_add_constraint(context.model, obj, MOI.ExponentialCone())
        end for i in 1:size(x, 1), j in 1:size(x, 2)
    ]
    context.constr_to_moi_inds[c] = inds
    return nothing
end

function populate_dual!(
    model::MOI.ModelLike,
    constr::ExpConstraint,
    MOI_constr_indices,
)
    x = first(children(constr))
    constr.dual = output([
        MOI.get(model, MOI.ConstraintDual(), MOI_constr_indices[i, j]) for
        i in 1:size(x, 1), j in 1:size(x, 2)
    ])
    return nothing
end
