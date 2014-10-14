import Base.isposdef
export SDPConstraint, isposdef

### Positive semidefinite cone constraint
type SDPConstraint <: Constraint
  head::Symbol
  children_hash::Uint64
  child::AbstractExpr
  size::(Int64, Int64)
  is_symmetric::Bool

  function SDPConstraint(child::AbstractExpr; is_symmetric=true)
    sz = child.size
    if sz[1] != sz[2]
      error("Positive semidefinite expressions must be square")
    end
    return new(:sdp, hash(child), child, sz, is_symmetric)
  end
end

function vexity(c::SDPConstraint)
  vex = vexity(c.child)
  if vex == AffineVexity() || vex == ConstVexity()
    return AffineVexity()
  else
    return NotDCP()
  end
end

function conic_form(c::SDPConstraint, unique_constr)
  if !((c.head, c.children_hash) in keys(unique_constr))
    objective, constraints = conic_form(c.child, unique_constr)
    new_constraint = ConicConstr([objective], :SDP, [c.size[1] * c.size[2]])
    push!(constraints, new_constraint)

    if c.is_symmetric
      _, new_constraints = conic_form(c.child == c.child', unique_constr)
      append!(constraints, new_constraints)
    end

    unique_constr[(c.head, c.children_hash)] = (objective, constraints)
  end
  return safe_copy(unique_constr[(c.head, c.children_hash)])
end

isposdef(x::AbstractExpr; is_symmetric=true) = SDPConstraint(x, is_symmetric=is_symmetric)

