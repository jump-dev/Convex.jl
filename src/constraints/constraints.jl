export EqConstraint, LtConstraint, GtConstraint
export ==, <=, >=

### Linear equality constraint
type EqConstraint <: Constraint
  head::Symbol
  id_hash::Uint64
  lhs::AbstractExpr
  rhs::AbstractExpr
  size::(Int, Int)

  function EqConstraint(lhs::AbstractExpr, rhs::AbstractExpr)
    if lhs.size == rhs.size || lhs.size == (1, 1)
      sz = rhs.size
    elseif rhs.size == (1, 1)
      sz = lhs.size
    else
      error("Cannot create equality constraint between expressions of size $(lhs.size) and $(rhs.size)")
    end
    return new(:(==), hash((lhs, rhs)), lhs, rhs, sz)
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

function conic_form!(c::EqConstraint, unique_conic_forms::UniqueConicForms)
  if !has_conic_form(unique_conic_forms, c)
    expr = c.lhs - c.rhs
    objective = conic_form!(expr, unique_conic_forms)
    new_constraint = ConicConstr([objective], :Zero, [c.size[1] * c.size[2]])
    cache_conic_form!(unique_conic_forms, c, new_constraint)
  end
  return get_conic_form(unique_conic_forms, c)
end

==(lhs::AbstractExpr, rhs::AbstractExpr) = EqConstraint(lhs, rhs)
==(lhs::AbstractExpr, rhs::Value) = ==(lhs, Constant(rhs))
==(lhs::Value, rhs::AbstractExpr) = ==(Constant(lhs), rhs)


### Linear inequality constraints
type LtConstraint <: Constraint
  head::Symbol
  id_hash::Uint64
  lhs::AbstractExpr
  rhs::AbstractExpr
  size::(Int, Int)

  function LtConstraint(lhs::AbstractExpr, rhs::AbstractExpr)
    if lhs.size == rhs.size || lhs.size == (1, 1)
      sz = rhs.size
    elseif rhs.size == (1, 1)
      sz = lhs.size
    else
      error("Cannot create inequality constraint between expressions of size $(lhs.size) and $(rhs.size)")
    end
    return new(:(<=), hash((lhs, rhs)), lhs, rhs, sz)
  end
end

function vexity(c::LtConstraint)
  vex = vexity(c.lhs) + (-vexity(c.rhs))
  if vex == ConcaveVexity()
    vex = NotDcp()
  end
  return vex
end

function conic_form!(c::LtConstraint, unique_conic_forms::UniqueConicForms)
  if !has_conic_form(unique_conic_forms, c)
    expr = c.rhs - c.lhs
    objective = conic_form!(expr, unique_conic_forms)
    new_constraint = ConicConstr([objective], :NonNeg, [c.size[1] * c.size[2]])
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


type GtConstraint <: Constraint
  head::Symbol
  id_hash::Uint64
  lhs::AbstractExpr
  rhs::AbstractExpr
  size::(Int, Int)

  function GtConstraint(lhs::AbstractExpr, rhs::AbstractExpr)
    if lhs.size == rhs.size || lhs.size == (1, 1)
      sz = rhs.size
    elseif rhs.size == (1, 1)
      sz = lhs.size
    else
      error("Cannot create inequality constraint between expressions of size $(lhs.size) and $(rhs.size)")
    end
    return new(:(>=), hash((lhs, rhs)), lhs, rhs, sz)
  end
end

function vexity(c::GtConstraint)
  vex = vexity(c.lhs) + (-vexity(c.rhs))
  if vex == ConvexVexity()
    vex = NotDcp()
  end
  return vex
end

function conic_form!(c::GtConstraint, unique_conic_forms::UniqueConicForms)
  if !has_conic_form(unique_conic_forms, c)
    expr = c.lhs - c.rhs
    objective = conic_form!(expr, unique_conic_forms)
    new_constraint = ConicConstr([objective], :NonNeg, [c.size[1] * c.size[2]])
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

function +{T<:Constraint, T2<:Constraint}(constraints_one::Array{T}, constraints_two::Array{T2})
  constraints = append!(Constraint[], constraints_one)
  return append!(constraints, constraints_two)
end
+(constraint_one::Constraint, constraint_two::Constraint) = [constraint_one] + [constraint_two]
+{T<:Constraint}(constraint_one::Constraint, constraints_two::Array{T}) =
  [constraint_one] + constraints_two
+{T<:Constraint}(constraints_one::Array{T}, constraint_two::Constraint) =
  constraints_one + [constraint_two]
