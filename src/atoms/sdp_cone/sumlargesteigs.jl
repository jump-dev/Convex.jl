#############################################################################
# sumlargesteigs.jl
# Handles maximum and minimum eigenvalue of a symmetric positive definite matrix
# (and imposes the constraint that its argument be PSD)
# All expressions and atoms are subtypes of AbstractExpr.
# Please read expressions.jl first.
#############################################################################
export sumlargesteigs

### Lambda max

type SumLargestEigs <: AbstractExpr
  head::Symbol
  id_hash::Uint64
  children::@compat Tuple{AbstractExpr}
  size::@compat Tuple{Int, Int}

  function SumLargestEigs(x::AbstractExpr, k::Int)
    children = (x, k)
    m,n = size(x)
    if m==n
      return new(:sumlargesteigs, hash(children), children, (1,1))
    else
      error("sumlargesteigs can only be applied to a square matrix.")
    end
  end
end

function sign(x::SumLargestEigs)
  return Positive()
end

function monotonicity(x::SumLargestEigs)
  return (Nondecreasing(),)
end

function curvature(x::SumLargestEigs)
  return ConvexVexity()
end

function evaluate(x::SumLargestEigs)
  eigvals(evaluate(x.children[1]))[end-x.children[2]:end]
end

sumlargesteigs(x::AbstractExpr) = SumLargestEigs(x)

# Create the equivalent conic problem ... ?  :
#   minimize sumlargest(l, k)
#   subject to
#            diagm(l) - A is positive semidefinite
#            A      is positive semidefinite
function conic_form!(x::SumLargestEigs, unique_conic_forms)
  if !has_conic_form(unique_conic_forms, x)
    A = x.children[1]
    k = x.children[2]
    m, n = size(A)
    l = Variable()
    p = minimize(sumlargest(l, k), diagm(l) - A ⪰ 0, A ⪰ 0)
    cache_conic_form!(unique_conic_forms, x, p)
  end
  return get_conic_form(unique_conic_forms, x)
end
