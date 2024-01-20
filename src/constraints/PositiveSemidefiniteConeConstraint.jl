mutable struct PositiveSemidefiniteConeConstraint <: Constraint
    child::AbstractExpr
    size::Tuple{Int,Int}
    dual::Union{Value,Nothing}

    function PositiveSemidefiniteConeConstraint(child::AbstractExpr)
        if child.size[1] != child.size[2]
            error("Positive semidefinite expressions must be square")
        end
        return new(child, child.size, nothing)
    end
end

head(io::IO, ::PositiveSemidefiniteConeConstraint) = print(io, "sdp")

AbstractTrees.children(c::PositiveSemidefiniteConeConstraint) = (c.child,)

function vexity(c::PositiveSemidefiniteConeConstraint)
    if !(vexity(c.child) in (AffineVexity(), ConstVexity()))
        return NotDcp()
    end
    return AffineVexity()
end

function _add_constraint!(
    context::Context,
    c::PositiveSemidefiniteConeConstraint,
)
    if vexity(c.child) == ConstVexity()
        x = evaluate(c.child)
        tol = CONSTANT_CONSTRAINT_TOL[]
        if !(x ≈ transpose(x))
            @warn "constant SDP constraint is violated"
            context.detected_infeasible_during_formulation[] = true
        elseif evaluate(LinearAlgebra.eigmin(c.child)) < -tol
            @warn "constant SDP constraint is violated"
            context.detected_infeasible_during_formulation[] = true
        end
        return
    end
    f = conic_form!(context, c.child)
    set = MOI.PositiveSemidefiniteConeSquare(c.size[1])
    context.constr_to_moi_inds[c] = MOI_add_constraint(context.model, f, set)
    return
end

function populate_dual!(
    model::MOI.ModelLike,
    c::PositiveSemidefiniteConeConstraint,
    indices,
)
    dual = MOI.get(model, MOI.ConstraintDual(), indices)
    c.dual = output(reshape(dual, c.size))
    return
end

function LinearAlgebra.isposdef(x::AbstractExpr)
    if iscomplex(x)
        return PositiveSemidefiniteConeConstraint(
            [real(x) -imag(x); imag(x) real(x)],
        )
    end
    return PositiveSemidefiniteConeConstraint(x)
end

⪰(x::AbstractExpr, y::AbstractExpr) = isposdef(x - y)

function ⪰(x::AbstractExpr, y::Value)
    if all(y .== 0)
        return isposdef(x)
    end
    return isposdef(x - constant(y))
end

function ⪰(x::Value, y::AbstractExpr)
    if all(x .== 0)
        return isposdef(-y)
    end
    return isposdef(constant(x) - y)
end

⪯(x::AbstractExpr, y::AbstractExpr) = ⪰(y, x)
⪯(x::Value, y::AbstractExpr) = ⪰(y, x)
⪯(x::AbstractExpr, y::Value) = ⪰(y, x)
