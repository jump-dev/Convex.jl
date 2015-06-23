#############################################################################
# sumlargesteigs.jl
# Handles top k eigenvalues of a symmetric positive definite matrix
# (and imposes the constraint that its argument be PSD)
# All expressions and atoms are subtypes of AbstractExpr.
# Please read expressions.jl first.
#############################################################################
export sumlargesteigs, schatten

### sumlargesteigs

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

# Create the equivalent conic problem:
#   minimize sk + Tr(Z)
#   subject to
#            Z ⪰ 0
#            A ⪰ 0
#            Z + sI ⪰ A
# See Ben-Tal and Nemirovski, "Lectures on Modern Convex Optimization"
# Example 18.c

function conic_form!(x::SumLargestEigs, unique_conic_forms)
  if !has_conic_form(unique_conic_forms, x)
    A = x.children[1]
    k = x.children[2]
    m, n = size(A)
    Z = Variable(n, n)
    s = Variable()
    p = minimize(s*k + trace(Z), 
                 Z + s*eye(n) - A ⪰ 0, 
                 A ⪰ 0, Z ⪰ 0)
    cache_conic_form!(unique_conic_forms, x, p)
  end
  return get_conic_form(unique_conic_forms, x)
end

### eigvals, for use in computing spectral functions

type EigAtom <: AbstractExpr
  head::Symbol
  id_hash::Uint64
  children::@compat Tuple{AbstractExpr}
  size::@compat Tuple{Int, Int}

  function EigAtom(x::AbstractExpr)
    children = (x,)
    m,n = size(x)
    if m==n
      return new(:eig, hash(children), children, (n,1))
    else
      error("eig can only be applied to a square matrix.")
    end
  end
end

function sign(x::EigAtom)
  return Positive()
end

# The monotonicity
function monotonicity(x::EigAtom)
  return (NoMonotonicity(),)
end

# eigvals can be composed with any symmetric function to give the right thing
# I guess we don't have a way to represent that...
function curvature(x::EigAtom)
  return AffineVexity()
end

function evaluate(x::EigAtom)
  return eigvals(evaluate(x.children[1]))
end

eigvals(x::AbstractExpr) = EigAtom(x)

# Create the equivalent conic constraints:
function conic_form!(x::EigAtom, unique_conic_forms)
  if !has_conic_form(unique_conic_forms, x)
    A = x.children[1]
    m, n = size(A)
    l = Variable(n)
    # monotonicity
    constraints = Constraint[l[i] >= l[i]+1 for i=1:n-1]
    # capture the right partial sums
    for k=1:n
      push!(constraints, sumlargesteigs(A,k) <= sum(l[1:k]))
    end  
    # and A had better be PSD
    push!(constraints, A ⪰ 0)
    conic_form!(constraints, unique_conic_forms)
    cache_conic_form!(unique_conic_forms, x, l)
  end
  return get_conic_form(unique_conic_forms, x)
end

schatten(x::AbstractExpr, p) = norm(eigenvalues(x), p)