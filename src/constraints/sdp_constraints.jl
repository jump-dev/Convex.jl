import Base.isposdef, Base.in
export SDPConstraint, isposdef, in

### Positive semidefinite cone constraint

# TODO: Terrible documentation. Please fix.
type SDPConstraint <: Constraint
  head::Symbol
  id_hash::Uint64
  child::AbstractExpr
  size::@compat Tuple{Int, Int}
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

function conic_form!(c::SDPConstraint, unique_conic_forms::UniqueConicForms)
  if !has_conic_form(unique_conic_forms, c)
    objective = conic_form!(c.child, unique_conic_forms)
    new_constraint = ConicConstr([objective], :SDP, [c.size[1] * c.size[2]])
    conic_constr_to_constr[new_constraint] = c
    cache_conic_form!(unique_conic_forms, c, new_constraint)
  end
  return get_conic_form(unique_conic_forms, c)
end

# TODO: Remove isposdef, change tests to use in. Update documentation and notebooks
isposdef(x::AbstractExpr) = SDPConstraint(x)

# TODO: Throw error if symbol is invalid.
function in(x::AbstractExpr, y::Symbol)
  if y == :semidefinite
    SDPConstraint(x)
  end
end
