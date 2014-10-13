import Base.isposdef
export SDPConstraint, isposdef

### Positive semidefinite cone constraint
type SDPConstraint <: Constraint
  head::Symbol
  children_hash::Uint64
  lhs::AbstractExpr
  size::(Int64, Int64)

  function SDPConstraint(lhs::AbstractExpr)
    sz = lhs.size
    if sz[1] != sz[2]
      error("Positive semidefinite expressions must be square")
    end
    return new(:sdp, hash(lhs), lhs, sz)
  end
end

function vexity(c::SDPConstraint)
  vex = vexity(c.lhs)
  if vex == AffineVexity() || vex == ConstVexity()
    return AffineVexity()
  else
    return NotDCP()
  end
end

function conic_form(c::SDPConstraint, unique_constr)
  if !((c.head, c.children_hash) in keys(unique_constr))
    objective, constraints = conic_form(c.lhs, unique_constr)
    new_constraint = ConicConstr([objective], :SDP, [c.size[1] * c.size[2]])
    push!(constraints, new_constraint)
    unique_constr[(c.head, c.children_hash)] = (objective, constraints)
  end
  return safe_copy(unique_constr[(c.head, c.children_hash)])
end

isposdef(x::AbstractExpr) = SDPConstraint(x)

