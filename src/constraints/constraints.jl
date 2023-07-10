import Base.==, Base.<=, Base.>=, Base.<, Base.>

# TODO- tighten this
const CONSTANT_CONSTRAINT_TOL = Ref(1e-2)

function iscomplex(constr::Constraint)
    return iscomplex(constr.lhs) || iscomplex(constr.rhs)
end

function add_constraint!(context::Context, c::Constraint)
    c ∈ keys(context.constr_to_moi_inds) && return
    return _add_constraint!(context, c)
end

### Linear equality constraint
mutable struct EqConstraint <: Constraint
    lhs::AbstractExpr
    rhs::AbstractExpr
    size::Tuple{Int,Int}
    dual::ValueOrNothing

    function EqConstraint(lhs::AbstractExpr, rhs::AbstractExpr)
        if lhs.size == rhs.size || lhs.size == (1, 1)
            sz = rhs.size
            if lhs.size == (1, 1) && rhs.size != (1, 1)
                lhs = lhs * ones(rhs.size)
            end
        elseif rhs.size == (1, 1)
            sz = lhs.size
            if rhs.size == (1, 1) && lhs.size != (1, 1)
                rhs = rhs * ones(lhs.size)
            end
        else
            error(
                "Cannot create equality constraint between expressions of size $(lhs.size) and $(rhs.size)",
            )
        end
        return new(lhs, rhs, sz, nothing)
    end
end

head(io::IO, ::EqConstraint) = print(io, "==")
function vexity(c::EqConstraint)
    vex = vexity(c.lhs) + (-vexity(c.rhs))
    # You can't have equality constraints with concave/convex expressions
    if vex == ConvexVexity() || vex == ConcaveVexity()
        vex = NotDcp()
    end
    return vex
end

function _add_constraint!(context::Context{T}, eq::EqConstraint) where {T}
    f = conic_form!(context, eq.lhs - eq.rhs)
    if f isa AbstractVector
        # a trivial constraint without variables like `5 == 0`
        if !all(abs.(f) .<= CONSTANT_CONSTRAINT_TOL[])
            @warn "Constant constraint is violated"
            context.detected_infeasible_during_formulation[] = true
        end
        return nothing
    end
    context.constr_to_moi_inds[eq] =
        MOI_add_constraint(context.model, f, MOI.Zeros(MOI.output_dimension(f)))
    return nothing
end

==(lhs::AbstractExpr, rhs::AbstractExpr) = EqConstraint(lhs, rhs)
==(lhs::AbstractExpr, rhs::Value) = ==(lhs, constant(rhs))
==(lhs::Value, rhs::AbstractExpr) = ==(constant(lhs), rhs)

### Linear inequality constraints
mutable struct LtConstraint <: Constraint
    lhs::AbstractExpr
    rhs::AbstractExpr
    size::Tuple{Int,Int}
    dual::ValueOrNothing

    function LtConstraint(lhs::AbstractExpr, rhs::AbstractExpr)
        if sign(lhs) == ComplexSign() || sign(rhs) == ComplexSign()
            error(
                "Cannot create inequality constraint between expressions of sign $(sign(lhs)) and $(sign(rhs))",
            )
        else
            if lhs.size == rhs.size || lhs.size == (1, 1)
                sz = rhs.size
                if lhs.size == (1, 1) && rhs.size != (1, 1)
                    lhs = lhs * ones(rhs.size)
                end
            elseif rhs.size == (1, 1)
                sz = lhs.size
                if rhs.size == (1, 1) && lhs.size != (1, 1)
                    rhs = rhs * ones(lhs.size)
                end
            else
                error(
                    "Cannot create inequality constraint between expressions of size $(lhs.size) and $(rhs.size)",
                )
            end
        end
        return new(lhs, rhs, sz, nothing)
    end
end
head(io::IO, ::LtConstraint) = print(io, "≤")

function vexity(c::LtConstraint)
    vex = vexity(c.lhs) + (-vexity(c.rhs))
    if vex == ConcaveVexity()
        vex = NotDcp()
    end
    return vex
end

function _add_constraint!(context::Context{T}, lt::LtConstraint) where {T}
    f = conic_form!(context, lt.rhs - lt.lhs)
    if f isa AbstractVector
        # a trivial constraint without variables like `5 >= 0`
        if !all(f .<= CONSTANT_CONSTRAINT_TOL[])
            @warn "Constant constraint is violated"
            context.detected_infeasible_during_formulation[] = true
        end
        return nothing
    end
    context.constr_to_moi_inds[lt] = MOI_add_constraint(
        context.model,
        f,
        MOI.Nonnegatives(MOI.output_dimension(f)),
    )
    return nothing
end

<=(lhs::AbstractExpr, rhs::AbstractExpr) = LtConstraint(lhs, rhs)
<=(lhs::AbstractExpr, rhs::Value) = <=(lhs, constant(rhs))
<=(lhs::Value, rhs::AbstractExpr) = <=(constant(lhs), rhs)
<(lhs::AbstractExpr, rhs::AbstractExpr) = LtConstraint(lhs, rhs)
<(lhs::AbstractExpr, rhs::Value) = <=(lhs, constant(rhs))
<(lhs::Value, rhs::AbstractExpr) = <=(constant(lhs), rhs)

mutable struct GtConstraint <: Constraint
    lhs::AbstractExpr
    rhs::AbstractExpr
    size::Tuple{Int,Int}
    dual::ValueOrNothing

    function GtConstraint(lhs::AbstractExpr, rhs::AbstractExpr)
        if sign(lhs) == ComplexSign() || sign(rhs) == ComplexSign()
            error(
                "Cannot create inequality constraint between expressions of sign $(sign(lhs)) and $(sign(rhs))",
            )
        else
            if lhs.size == rhs.size || lhs.size == (1, 1)
                sz = rhs.size
                if lhs.size == (1, 1) && rhs.size != (1, 1)
                    lhs = lhs * ones(rhs.size)
                end
            elseif rhs.size == (1, 1)
                sz = lhs.size
                if rhs.size == (1, 1) && lhs.size != (1, 1)
                    rhs = rhs * ones(lhs.size)
                end
            else
                error(
                    "Cannot create inequality constraint between expressions of size $(lhs.size) and $(rhs.size)",
                )
            end
        end
        return new(lhs, rhs, sz, nothing)
    end
end
head(io::IO, ::GtConstraint) = print(io, "≥")

function vexity(c::GtConstraint)
    vex = -vexity(c.lhs) + (vexity(c.rhs))
    if vex == ConcaveVexity()
        vex = NotDcp()
    end
    return vex
end

function _add_constraint!(context::Context{T}, gt::GtConstraint) where {T}
    f = conic_form!(context, gt.lhs - gt.rhs)
    if f isa AbstractVector
        # a trivial constraint without variables like `5 >= 0`
        if !all(f .>= -CONSTANT_CONSTRAINT_TOL[])
            @warn "Constant constraint is violated"
            context.detected_infeasible_during_formulation[] = true
        end
        return nothing
    end
    context.constr_to_moi_inds[gt] = MOI_add_constraint(
        context.model,
        f,
        MOI.Nonnegatives(MOI.output_dimension(f)),
    )
    return nothing
end

>=(lhs::AbstractExpr, rhs::AbstractExpr) = GtConstraint(lhs, rhs)
>=(lhs::AbstractExpr, rhs::Value) = >=(lhs, constant(rhs))
>=(lhs::Value, rhs::AbstractExpr) = >=(constant(lhs), rhs)
>(lhs::AbstractExpr, rhs::AbstractExpr) = GtConstraint(lhs, rhs)
>(lhs::AbstractExpr, rhs::Value) = >=(lhs, constant(rhs))
>(lhs::Value, rhs::AbstractExpr) = >=(constant(lhs), rhs)

function +(
    constraints_one::Array{<:Constraint},
    constraints_two::Array{<:Constraint},
)
    constraints = append!(Constraint[], constraints_one)
    return append!(constraints, constraints_two)
end
function +(constraint_one::Constraint, constraint_two::Constraint)
    return [constraint_one] + [constraint_two]
end
function +(constraint_one::Constraint, constraints_two::Array{<:Constraint})
    return [constraint_one] + constraints_two
end
function +(constraints_one::Array{<:Constraint}, constraint_two::Constraint)
    return constraints_one + [constraint_two]
end

function populate_dual!(
    model::MOI.ModelLike,
    constr::Union{EqConstraint,GtConstraint,LtConstraint},
    MOI_constr_indices,
)
    if iscomplex(constr)
        re = MOI.get(model, MOI.ConstraintDual(), MOI_constr_indices[1])
        imag = MOI.get(model, MOI.ConstraintDual(), MOI_constr_indices[2])
        constr.dual = output(reshape(re + im * imag, constr.size))
    else
        constr.dual = output(
            reshape(
                MOI.get(model, MOI.ConstraintDual(), MOI_constr_indices),
                constr.size,
            ),
        )
    end
end
