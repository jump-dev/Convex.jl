import Base.isposdef
export SDPConstraint, isposdef

### Positive semidefinite cone constraint
type SDPConstraint <: Constraint
  head::Symbol
  id_hash::Uint64
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

function conic_form!(c::SDPConstraint, unique_conic_forms::UniqueConicForms)
  if !has_conic_form(unique_conic_forms, c)
    objective = conic_form!(c.child, unique_conic_forms)
    if c.is_symmetric
      n,m = size(c.child)
      for i=1:n
        for j=i+1:m
          # it would be worthwhile to make this more efficient by vectorizing
          conic_form!(c.child[i,j] - c.child[j,i]==0, unique_conic_forms)
        end
      end
    end
    cache_conic_form!(unique_conic_forms, c, ConicConstr([objective], :SDP, [c.size[1] * c.size[2]]))
  end
  return get_conic_form(unique_conic_forms, c)
end

isposdef(x::AbstractExpr; is_symmetric=true) = SDPConstraint(x, is_symmetric=is_symmetric)

