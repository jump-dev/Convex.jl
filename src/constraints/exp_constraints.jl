### (Primal) exponential cone constraint ExpConstraint(x,y,z) => y exp(x/y) <= z & y>=0
struct ExpConstraint <: Constraint
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
    return ExpConstraint(x, Constant(y), z)
end
function ExpConstraint(x::AbstractExpr, y::AbstractExpr, z)
    return ExpConstraint(x, y, Constant(z))
end
function ExpConstraint(x, y::AbstractExpr, z::AbstractExpr)
    return ExpConstraint(Constant(x), y, z)
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

function conic_form!(c::ExpConstraint, unique_conic_forms::UniqueConicForms)
    if !has_conic_form(unique_conic_forms, c)
        conic_constrs = ConicConstr[]
        if c.size == (1, 1)
            objectives = Vector{ConicObj}(undef, 3)
            @inbounds for iobj in 1:3
                objectives[iobj] =
                    conic_form!(c.children[iobj], unique_conic_forms)
            end
            push!(conic_constrs, ConicConstr(objectives, :ExpPrimal, [1, 1, 1]))
        else
            for i in 1:c.size[1]
                for j in 1:c.size[2]
                    objectives = Vector{ConicObj}(undef, 3)
                    @inbounds for iobj in 1:3
                        objectives[iobj] = conic_form!(
                            c.children[iobj][i, j],
                            unique_conic_forms,
                        )
                    end
                    push!(
                        conic_constrs,
                        ConicConstr(objectives, :ExpPrimal, [1, 1, 1]),
                    )
                end
            end
        end
        cache_conic_form!(unique_conic_forms, c, conic_constrs)
    end
    return get_conic_form(unique_conic_forms, c)
end
