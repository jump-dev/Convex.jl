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

# users specify SDPs as `A in :SDP` where A is an n x n square matrix
# solvers (Mosek and SCS) specify the *lower triangular part of A* is in the SDP cone
# so we need the lower triangular part (n*(n+1)/2 entries) of c.child to be in the SDP cone
# and we need the corresponding upper elements to match the lower elements, 
# which we enforce via equality constraints
function conic_form!(c::SDPConstraint, unique_conic_forms::UniqueConicForms)
  if !has_conic_form(unique_conic_forms, c)
    # construct linear indices to pick out the lower triangular part (including diagonal),
    # the upper triangular part (not including diagonal)
    # and the corresponding entries in the lower triangular part, so
    # symmetry => c.child[upperpart] 
    lowerpart = Array(Int, int(n*(n-1)/2))
    upperpart = Array(Int, int(n*(n-1)/2))
    diagandlowerpart = Array(Int, int(n*(n+1)/2))
    kdiag, klower = 0, 0
    for i = 1:n
      for j = 1:i
        if j < i # on the strictly lower part
          klower += 1
          diagandlowerpart[kdiag + klower] = n*(j-1) + i
          upperpart[klower] = n*(i-1) + j
          lowerpart[klower] = n*(j-1) + i
        else # on the diagonal
          kdiag += 1 
          diagandlowerpart[kdiag + klower] = n*(j-1) + i
        end
      end
    end
    objective = conic_form!(c.child[diagandlowerpart], unique_conic_forms)
    new_constraint = ConicConstr([objective], :SDP, [int(n*(n+1)/2)])
    conic_constr_to_constr[new_constraint] = c
    cache_conic_form!(unique_conic_forms, c, new_constraint)

    # make sure upper and lower triangular part match in the solution
    # TODO: 1) propagate dual values 2) presolve to eliminate these variables
    equality_constraint = conic_form!(c.child[lowerpart] == c.child[upperpart], unique_conic_forms)
  end
  return get_conic_form(unique_conic_forms, c)
end

# TODO: Remove isposdef, change tests to use is. Update documentation and notebooks
isposdef(x::AbstractExpr) = SDPConstraint(x)

# TODO: Throw error if symbol is invalid.
function in(x::AbstractExpr, y::Symbol)
  if y == :semidefinite
    SDPConstraint(x)
  end
end
