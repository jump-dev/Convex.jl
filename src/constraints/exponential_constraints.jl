### (Primal) exponential cone constraint ExpConstraint(x,y,z) => y exp(x/y) <= z
type ExpConstraint <: Constraint
  head::Symbol
  children_hash::Uint64
  objective::Integer # index of the child that represents the objective
  children::(AbstractExpr, AbstractExpr, AbstractExpr) # (x, y, z)
  size::(Int64, Int64)

  function ExpConstraint(objective::Integer, x::AbstractExpr, y::AbstractExpr, z::AbstractExpr)
    @assert(x.size == y.size == z.size,
           "Exponential constraint requires x, y, and z to be of same size")
    # @assert(x.size == (1,1),
    #        "Exponential constraint requires x, y, and z to be scalar for now")
    sz = x.size
    return new(:exp, hash((x,y,z)), objective, (x, y, z), sz)
  end
end
ExpConstraint(objective::Integer, x::AbstractExpr, y, z::AbstractExpr) = ExpConstraint(objective, x, Constant(y), z)

function vexity(c::ExpConstraint)
  # TODO: check these...
  if vexity(c.x) == ConcaveVexity()
    error("Exponential constraint requires x to be convex")
  end
  if vexity(c.y)!=ConstVexity()
    error("Exponential constraint requires y to be constant")
  end
  if vexity(c.z) == ConvexVexity()
    error("Exponential constraint requires z to be concave")
  end
  return ConvexVexity()
end

function conic_form(c::ExpConstraint, unique_constr)
  if !((c.head, c.children_hash) in keys(unique_constr))
    constraints = ConicConstr[]
    objectives = Array(ConicObj, 3)
    for iobj=1:3
      objectives[iobj], new_constraints = conic_form(c.children[iobj], unique_constr)
      append!(constraints, new_constraints)
    end
    unique_constr[(c.head, c.children_hash)] = (objectives[c.objective], constraints)
    if c.size == (1, 1)
      exp_constraint = ConicConstr(objectives, :ExpPrimal, [1, 1, 1])
      push!(constraints, exp_constraint)
    else
      for i=1:c.size[1]
        for j=1:c.size[2]
          objectives = Array(ConicObj, 3)
          for iobj=1:3
            objectives[iobj], _ = conic_form(c.children[iobj][i,j], unique_constr)
          end
          exp_constraint = ConicConstr(objectives, :ExpPrimal, [1, 1, 1])
          push!(constraints, exp_constraint)
        end
      end
    end
  end
  return safe_copy(unique_constr[(c.head, c.children_hash)])
end
