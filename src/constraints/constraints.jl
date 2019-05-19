import Base.==, Base.<=, Base.>=, Base.<, Base.>
export EqConstraint, LtConstraint, GtConstraint
export ==, <=, >=

const conic_constr_to_constr = Dict{ConicConstr, Constraint}()

### Linear equality constraint
mutable struct EqConstraint <: Constraint
    head::Symbol
    id_hash::UInt64
    lhs::AbstractExpr
    rhs::AbstractExpr
    size::Tuple{Int, Int}
    dual::ValueOrNothing

    function EqConstraint(lhs::AbstractExpr, rhs::AbstractExpr)
        if lhs.size == rhs.size || lhs.size == (1, 1)
            sz = rhs.size
        elseif rhs.size == (1, 1)
            sz = lhs.size
        else
            error("cannot create equality constraint between expressions of size $(lhs.size) and $(rhs.size)")
        end
        id_hash = hash((lhs, rhs, :(==)))
        return new(:(==), id_hash, lhs, rhs, sz, nothing)
    end
end

function vexity(c::EqConstraint)
    vex = vexity(c.lhs) + (-vexity(c.rhs))
    # You can't have equality constraints with concave/convex expressions
    if vex == ConvexVexity() || vex == ConcaveVexity()
        vex = NotDcp()
    end
    return vex
end

function conic_form!(c::EqConstraint, unique_conic_forms::UniqueConicForms=UniqueConicForms())
    if !has_conic_form(unique_conic_forms, c)
        if !(sign(c.lhs) == ComplexSign() || sign(c.rhs) == ComplexSign())

            expr = c.lhs - c.rhs
            objective = conic_form!(expr, unique_conic_forms)
            new_constraint = ConicConstr([objective], :Zero, [c.size[1] * c.size[2]])
            conic_constr_to_constr[new_constraint] = c
        else
            real_expr = real(c.lhs - c.rhs)
            imag_expr = imag(c.lhs - c.rhs)
            real_objective = conic_form!(real_expr, unique_conic_forms)
            imag_objective = conic_form!(imag_expr, unique_conic_forms)
            new_constraint = ConicConstr([real_objective, imag_objective], :Zero, [c.size[1] * c.size[2], c.size[1] * c.size[2]])
            conic_constr_to_constr[new_constraint] = c
        end
        cache_conic_form!(unique_conic_forms, c, new_constraint)
    end
    return get_conic_form(unique_conic_forms, c)
end

==(lhs::AbstractExpr, rhs::AbstractExpr) = EqConstraint(lhs, rhs)
==(lhs::AbstractExpr, rhs::Value) = ==(lhs, Constant(rhs))
==(lhs::Value, rhs::AbstractExpr) = ==(Constant(lhs), rhs)


### Linear inequality constraints
mutable struct LtConstraint <: Constraint
    head::Symbol
    id_hash::UInt64
    lhs::AbstractExpr
    rhs::AbstractExpr
    size::Tuple{Int, Int}
    dual::ValueOrNothing

    function LtConstraint(lhs::AbstractExpr, rhs::AbstractExpr)
        if sign(lhs) == ComplexSign() || sign(rhs) == ComplexSign()
            error("cannot create inequality constraint between expressions of sign $(sign(lhs)) and $(sign(rhs))")
        else
            if lhs.size == rhs.size || lhs.size == (1, 1)
                sz = rhs.size
            elseif rhs.size == (1, 1)
                sz = lhs.size
            else
                error("cannot create inequality constraint between expressions of size $(lhs.size) and $(rhs.size)")
            end
        end
        id_hash = hash((lhs, rhs, :(<=)))
        return new(:(<=), id_hash, lhs, rhs, sz, nothing)
    end
end

function vexity(c::LtConstraint)
    vex = vexity(c.lhs) + (-vexity(c.rhs))
    if vex == ConcaveVexity()
        vex = NotDcp()
    end
    return vex
end

function conic_form!(c::LtConstraint, unique_conic_forms::UniqueConicForms=UniqueConicForms())
    if !has_conic_form(unique_conic_forms, c)
        expr = c.rhs - c.lhs
        objective = conic_form!(expr, unique_conic_forms)
        new_constraint = ConicConstr([objective], :NonNeg, [c.size[1] * c.size[2]])
        conic_constr_to_constr[new_constraint] = c
        cache_conic_form!(unique_conic_forms, c, new_constraint)
    end
    return get_conic_form(unique_conic_forms, c)
end

<=(lhs::AbstractExpr, rhs::AbstractExpr) = LtConstraint(lhs, rhs)
<=(lhs::AbstractExpr, rhs::Value) = <=(lhs, Constant(rhs))
<=(lhs::Value, rhs::AbstractExpr) = <=(Constant(lhs), rhs)
<(lhs::AbstractExpr, rhs::AbstractExpr) = LtConstraint(lhs, rhs)
<(lhs::AbstractExpr, rhs::Value) = <=(lhs, Constant(rhs))
<(lhs::Value, rhs::AbstractExpr) = <=(Constant(lhs), rhs)


mutable struct GtConstraint <: Constraint
    head::Symbol
    id_hash::UInt64
    lhs::AbstractExpr
    rhs::AbstractExpr
    size::Tuple{Int, Int}
    dual::ValueOrNothing

    function GtConstraint(lhs::AbstractExpr, rhs::AbstractExpr)
        if sign(lhs) == ComplexSign() || sign(rhs) == ComplexSign()
            error("cannot create inequality constraint between expressions of sign $(sign(lhs)) and $(sign(rhs))")
        else
            if lhs.size == rhs.size || lhs.size == (1, 1)
                sz = rhs.size
            elseif rhs.size == (1, 1)
                sz = lhs.size
            else
                error("cannot create inequality constraint between expressions of size $(lhs.size) and $(rhs.size)")
            end
        end
        id_hash = hash((lhs, rhs, :(>=)))
        return new(:(>=), id_hash, lhs, rhs, sz, nothing)
    end
end

function vexity(c::GtConstraint)
    vex = -vexity(c.lhs) + (vexity(c.rhs))
    if vex == ConcaveVexity()
        vex = NotDcp()
    end
    return vex
end

function conic_form!(c::GtConstraint, unique_conic_forms::UniqueConicForms=UniqueConicForms())
    if !has_conic_form(unique_conic_forms, c)
        expr = c.lhs - c.rhs
        objective = conic_form!(expr, unique_conic_forms)
        new_constraint = ConicConstr([objective], :NonNeg, [c.size[1] * c.size[2]])
        conic_constr_to_constr[new_constraint] = c
        cache_conic_form!(unique_conic_forms, c, new_constraint)
    end
    return get_conic_form(unique_conic_forms, c)
end

>=(lhs::AbstractExpr, rhs::AbstractExpr) = GtConstraint(lhs, rhs)
>=(lhs::AbstractExpr, rhs::Value) = >=(lhs, Constant(rhs))
>=(lhs::Value, rhs::AbstractExpr) = >=(Constant(lhs), rhs)
>(lhs::AbstractExpr, rhs::AbstractExpr) = GtConstraint(lhs, rhs)
>(lhs::AbstractExpr, rhs::Value) = >=(lhs, Constant(rhs))
>(lhs::Value, rhs::AbstractExpr) = >=(Constant(lhs), rhs)

function +(constraints_one::Array{<:Constraint}, constraints_two::Array{<:Constraint})
    constraints = append!(Constraint[], constraints_one)
    return append!(constraints, constraints_two)
end
+(constraint_one::Constraint, constraint_two::Constraint) = [constraint_one] + [constraint_two]
+(constraint_one::Constraint, constraints_two::Array{<:Constraint}) =
    [constraint_one] + constraints_two
+(constraints_one::Array{<:Constraint}, constraint_two::Constraint) =
    constraints_one + [constraint_two]
