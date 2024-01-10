mutable struct SDPConstraint <: Constraint
    child::AbstractExpr
    size::Tuple{Int,Int}
    dual::Union{Value,Nothing}

    function SDPConstraint(child::AbstractExpr)
        if child.size[1] != child.size[2]
            error("Positive semidefinite expressions must be square")
        end
        return new(child, child.size, nothing)
    end
end

head(io::IO, ::SDPConstraint) = print(io, "sdp")

function vexity(c::SDPConstraint)
    vex = vexity(c.child)
    if vex == AffineVexity() || vex == ConstVexity()
        return AffineVexity()
    end
    return NotDcp()
end

function _add_constraint!(context::Context, c::SDPConstraint)
    if vexity(c.child) == ConstVexity()
        x = evaluate(c.child)
        if !(x ≈ transpose(x))
            @warn "constant SDP constraint is violated"
            context.detected_infeasible_during_formulation[] = true
        end
        if evaluate(LinearAlgebra.eigmin(c.child)) < -CONSTANT_CONSTRAINT_TOL[]
            @warn "constant SDP constraint is violated"
            context.detected_infeasible_during_formulation[] = true
        end
        return
    end
    context.constr_to_moi_inds[c] = MOI_add_constraint(
        context.model,
        conic_form!(context, c.child),
        MOI.PositiveSemidefiniteConeSquare(c.size[1]),
    )
    return
end

function populate_dual!(model::MOI.ModelLike, c::SDPConstraint, indices)
    dual = MOI.get(model, MOI.ConstraintDual(), indices)
    c.dual = output(reshape(dual, c.size))
    return
end

# TODO: Remove isposdef, change tests to use in. Update documentation and
# notebooks
LinearAlgebra.isposdef(x::AbstractExpr) = in(x, :SDP)

function Base.in(x::AbstractExpr, y::Symbol)
    if !(y in (:semidefinite, :SDP))
        error("Set $y not understood")
    end
    if iscomplex(x)
        return SDPConstraint([real(x) -imag(x); imag(x) real(x)])
    end
    return SDPConstraint(x)
end

⪰(x::AbstractExpr, y::AbstractExpr) = in(x - y, :SDP)

function ⪰(x::AbstractExpr, y::Value)
    if all(y .== 0)
        return in(x, :SDP)
    end
    return in(x - constant(y), :SDP)
end

function ⪰(x::Value, y::AbstractExpr)
    if all(x .== 0)
        return in(-y, :SDP)
    end
    return in(constant(x) - y, :SDP)
end

⪯(x::AbstractExpr, y::AbstractExpr) = ⪰(y, x)
⪯(x::Value, y::AbstractExpr) = ⪰(y, x)
⪯(x::AbstractExpr, y::Value) = ⪰(y, x)
