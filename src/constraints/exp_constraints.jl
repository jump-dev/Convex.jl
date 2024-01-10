mutable struct ExpConstraint <: Constraint
    children::Tuple{AbstractExpr,AbstractExpr,AbstractExpr} # (x, y, z)
    size::Tuple{Int,Int}
    dual::Union{Value,Nothing}

    function ExpConstraint(x::AbstractExpr, y::AbstractExpr, z::AbstractExpr)
        @assert(
            x.size == y.size == z.size,
            "Exponential constraint requires x, y, and z to be of same size"
        )
        return new((x, y, z), x.size, nothing)
    end
end

head(io::IO, ::ExpConstraint) = print(io, "exp")

function ExpConstraint(x::AbstractExpr, y, z::AbstractExpr)
    return ExpConstraint(x, constant(y), z)
end

function ExpConstraint(x::AbstractExpr, y::AbstractExpr, z)
    return ExpConstraint(x, y, constant(z))
end

function ExpConstraint(x, y::AbstractExpr, z::AbstractExpr)
    return ExpConstraint(constant(x), y, z)
end

# function vexity(c::ExpConstraint)
#     # TODO: check these...
#     if vexity(c.children[1]) == ConcaveVexity()
#         error("Exponential constraint requires x to be convex")
#     end
#     if vexity(c.children[2]) != ConstVexity()
#         error("Exponential constraint requires y to be constant")
#     end
#     if vexity(c.children[3]) == ConvexVexity()
#         error("Exponential constraint requires z to be concave")
#     end
#     return ConvexVexity()
# end

function _add_constraint!(context::Context{T}, c::ExpConstraint) where {T}
    context.constr_to_moi_inds[c] = [
        MOI_add_constraint(
            context.model,
            operate(
                vcat,
                T,
                sum(map(sign, c.children)),
                conic_form!(context, c.children[1][i, j]),
                conic_form!(context, c.children[2][i, j]),
                conic_form!(context, c.children[3][i, j]),
            ),
            MOI.ExponentialCone(),
        ) for i in 1:size(c.children[1], 1), j in 1:size(c.children[1], 2)
    ]
    return
end

function populate_dual!(model::MOI.ModelLike, c::ExpConstraint, indices)
    c.dual = output(MOI.get.(model, MOI.ConstraintDual(), indices))
    return
end
