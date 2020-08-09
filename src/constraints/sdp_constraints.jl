import LinearAlgebra.isposdef
import Base.in
### Positive semidefinite cone constraint

# TODO: Terrible documentation. Please fix.
mutable struct SDPConstraint <: Constraint
    head::Symbol
    id_hash::UInt64
    child::AbstractExpr
    size::Tuple{Int, Int}
    dual::ValueOrNothing

    function SDPConstraint(child::AbstractExpr)
        sz = child.size
        if sz[1] != sz[2]
            error("Positive semidefinite expressions must be square")
        end
        id_hash = hash((child, :sdp))
        return new(:sdp, id_hash, child, sz, nothing)
    end
end

function vexity(c::SDPConstraint)
    vex = vexity(c.child)
    if vex == AffineVexity() || vex == ConstVexity()
        return AffineVexity()
    else
        return NotDcp()
    end
end

function _add_constraints_to_context(c::SDPConstraint, context::Context)
    f = template(c.child, context)
    d = c.size[1]
    context.constr_to_moi_inds[c] = MOI_add_constraint(context.model, f, MOI.PositiveSemidefiniteConeSquare(d))
    return nothing
end

function populate_dual!(model::MOI.ModelLike, constr::SDPConstraint, MOI_constr_indices)
    constr.dual = output(reshape(MOI.get(model, MOI.ConstraintDual(), MOI_constr_indices), constr.size))
end


# TODO: Remove isposdef, change tests to use in. Update documentation and notebooks
function isposdef(x::AbstractExpr)
    if iscomplex(x)
        SDPConstraint([real(x) -imag(x);imag(x) real(x)])
    else
        SDPConstraint(x)
    end
end

function in(x::AbstractExpr, y::Symbol)
    if y == :semidefinite || y == :SDP
        if iscomplex(x)
            SDPConstraint([real(x) -imag(x);imag(x) real(x)])
        else
            SDPConstraint(x)
        end
    else
        error("Set $y not understood")
    end
end

function ⪰(x::AbstractExpr, y::AbstractExpr)
    if iscomplex(x) || iscomplex(y)
        SDPConstraint([real(x-y) -imag(x-y);imag(x-y) real(x-y)])
    else
        SDPConstraint(x-y)
    end
end

function ⪯(x::AbstractExpr, y::AbstractExpr)
    if iscomplex(x) || iscomplex(y)
        SDPConstraint([real(y-x) -imag(y-x);imag(y-x) real(y-x)])
    else
        SDPConstraint(y-x)
    end
end

function ⪰(x::AbstractExpr, y::Value)
    if iscomplex(x) || iscomplex(y)
        all(y .== 0) ? SDPConstraint([real(x) -imag(x);imag(x) real(x)]) : SDPConstraint([real(x-constant(y)) -imag(x-constant(y));imag(x-constant(y)) real(x-constant(y))])
    else
        all(y .== 0) ? SDPConstraint(x) : SDPConstraint(x - constant(y))
    end

end

function ⪰(x::Value, y::AbstractExpr)
    if iscomplex(y) || iscomplex(x)
        all(x .== 0) ? SDPConstraint([real(-y) -imag(-y);imag(-y) real(-y)]) : SDPConstraint([real(constant(x)-y) -imag(constant(x)-y);imag(constant(x)-y) real(constant(x)-y)])
    else
        all(x .== 0) ? SDPConstraint(-y) : SDPConstraint(constant(x) - y)
    end
end

function ⪯(x::Value, y::AbstractExpr)
    if iscomplex(y) || iscomplex(x)
        all(x .== 0) ? SDPConstraint([real(y) -imag(y);imag(y) real(y)]) : SDPConstraint([real(y-constant(x)) -imag(y-constant(x));imag(y-constant(x)) real(y-constant(x))])
    else
        all(x .== 0) ? SDPConstraint(y) : SDPConstraint(y - constant(x))
    end
end

function ⪯(x::AbstractExpr, y::Value)
    if iscomplex(x) || iscomplex(y)
        all(y .== 0) ? SDPConstraint([real(-x) -imag(-x);imag(-x) real(-x)]) : SDPConstraint([real(constant(y)-x) -imag(constant(y)-x);imag(constant(y)-x) real(constant(y)-x)])
    else
        all(y .== 0) ? SDPConstraint(-x) : SDPConstraint(constant(y) - x)
    end
end
