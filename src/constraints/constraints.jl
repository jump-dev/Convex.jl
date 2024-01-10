const CONSTANT_CONSTRAINT_TOL = Ref(1e-6)

Base.:+(x::Array{<:Constraint}, y::Array{<:Constraint}) = vcat(x, y)

Base.:+(x::Constraint, y::Constraint) = [x, y]

Base.:+(x::Constraint, y::Array{<:Constraint}) = vcat(x, y)

Base.:+(x::Array{<:Constraint}, y::Constraint) = vcat(x, y)

iscomplex(c::Constraint) = iscomplex(c.lhs) || iscomplex(c.rhs)

function add_constraint!(context::Context, c::Constraint)
    if c in keys(context.constr_to_moi_inds)
        return
    end
    return _add_constraint!(context, c)
end

mutable struct EqConstraint <: Constraint
    lhs::AbstractExpr
    rhs::AbstractExpr
    size::Tuple{Int,Int}
    dual::Union{Value,Nothing}

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
        return NotDcp()
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
        return
    end
    context.constr_to_moi_inds[eq] =
        MOI_add_constraint(context.model, f, MOI.Zeros(MOI.output_dimension(f)))
    return
end

Base.:(==)(lhs::AbstractExpr, rhs::AbstractExpr) = EqConstraint(lhs, rhs)

Base.:(==)(lhs::AbstractExpr, rhs::Value) = ==(lhs, constant(rhs))

Base.:(==)(lhs::Value, rhs::AbstractExpr) = ==(constant(lhs), rhs)

mutable struct LtConstraint <: Constraint
    lhs::AbstractExpr
    rhs::AbstractExpr
    size::Tuple{Int,Int}
    dual::Union{Value,Nothing}

    function LtConstraint(lhs::AbstractExpr, rhs::AbstractExpr)
        if sign(lhs) == ComplexSign() || sign(rhs) == ComplexSign()
            error(
                "Cannot create inequality constraint between expressions of sign $(sign(lhs)) and $(sign(rhs))",
            )
        end
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
        return new(lhs, rhs, sz, nothing)
    end
end

head(io::IO, ::LtConstraint) = print(io, "≤")

function vexity(c::LtConstraint)
    vex = vexity(c.lhs) + (-vexity(c.rhs))
    if vex == ConcaveVexity()
        return NotDcp()
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
        return
    end
    context.constr_to_moi_inds[lt] = MOI_add_constraint(
        context.model,
        f,
        MOI.Nonnegatives(MOI.output_dimension(f)),
    )
    return
end

Base.:<=(lhs::AbstractExpr, rhs::AbstractExpr) = LtConstraint(lhs, rhs)

Base.:<=(lhs::AbstractExpr, rhs::Value) = <=(lhs, constant(rhs))

Base.:<=(lhs::Value, rhs::AbstractExpr) = <=(constant(lhs), rhs)

Base.:<(lhs::AbstractExpr, rhs::AbstractExpr) = LtConstraint(lhs, rhs)

Base.:<(lhs::AbstractExpr, rhs::Value) = <=(lhs, constant(rhs))

Base.:<(lhs::Value, rhs::AbstractExpr) = <=(constant(lhs), rhs)

mutable struct GtConstraint <: Constraint
    lhs::AbstractExpr
    rhs::AbstractExpr
    size::Tuple{Int,Int}
    dual::Union{Value,Nothing}

    function GtConstraint(lhs::AbstractExpr, rhs::AbstractExpr)
        if sign(lhs) == ComplexSign() || sign(rhs) == ComplexSign()
            error(
                "Cannot create inequality constraint between expressions of sign $(sign(lhs)) and $(sign(rhs))",
            )
        end
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
        return new(lhs, rhs, sz, nothing)
    end
end

head(io::IO, ::GtConstraint) = print(io, "≥")

function vexity(c::GtConstraint)
    vex = -vexity(c.lhs) + (vexity(c.rhs))
    if vex == ConcaveVexity()
        return NotDcp()
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
        return
    end
    context.constr_to_moi_inds[gt] = MOI_add_constraint(
        context.model,
        f,
        MOI.Nonnegatives(MOI.output_dimension(f)),
    )
    return
end

Base.:>=(lhs::AbstractExpr, rhs::AbstractExpr) = GtConstraint(lhs, rhs)

Base.:>=(lhs::AbstractExpr, rhs::Value) = >=(lhs, constant(rhs))

Base.:>=(lhs::Value, rhs::AbstractExpr) = >=(constant(lhs), rhs)

Base.:>(lhs::AbstractExpr, rhs::AbstractExpr) = GtConstraint(lhs, rhs)

Base.:>(lhs::AbstractExpr, rhs::Value) = >=(lhs, constant(rhs))

Base.:>(lhs::Value, rhs::AbstractExpr) = >=(constant(lhs), rhs)

function populate_dual!(
    model::MOI.ModelLike,
    c::Union{EqConstraint,GtConstraint,LtConstraint},
    indices,
)
    if iscomplex(c)
        re = MOI.get(model, MOI.ConstraintDual(), indices[1])
        imag = MOI.get(model, MOI.ConstraintDual(), indices[2])
        c.dual = output(reshape(re + im * imag, c.size))
    else
        ret = MOI.get(model, MOI.ConstraintDual(), indices)
        c.dual = output(reshape(ret, c.size))
    end
    return
end
