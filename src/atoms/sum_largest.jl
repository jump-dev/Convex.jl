#############################################################################
# sum_largest.jl
# The sum of the k largest/smallest elements of an expression
# All expressions and atoms are subtpyes of AbstractExpr.
# Please read expressions.jl first.
#############################################################################

export sum_largest, sum_smallest, weighted_sum
export sign, curvature, monotonicity, evaluate

type SumLargestAtom <: AbstractExpr
  head::Symbol
  id_hash::Uint64
  children::(AbstractExpr,)
  size::(Int, Int)
  k::Int

  function SumLargestAtom(x::AbstractExpr, k::Int)
    if k <= 0
      error("sum_largest and sum_smallest only support positive values of k")
    end
    if k > get_vectorized_size(x)
      error("k cannot be larger than the number of entries in x")
    end
    children = (x,)
    return new(:sum_largest, hash((children, k)), children, (1,1), k)
  end
end

function sign(x::SumLargestAtom)
  return sign(x.children[1])
end

function monotonicity(x::SumLargestAtom)
  return (Nondecreasing(), )
end

function curvature(x::SumLargestAtom)
  return ConvexVexity()
end

function evaluate(x::SumLargestAtom)
  return sum(sort(vec(evaluate(x.children[1])), rev=true)[1:x.k])
end

function conic_form!(x::SumLargestAtom, unique_conic_forms::UniqueConicForms)
  if !has_conic_form(unique_conic_forms, x)
    c = x.children[1]
    t = Variable(size(c))
    q = Variable()
    # sum k largest given by the solution to
    # minimize sum(t) + k*q
    # subject to c <= t + q, t >= 0
    objective = conic_form!(sum(t) + x.k*q, unique_conic_forms)
    conic_form!(c <= t + q, unique_conic_forms)
    conic_form!(t >= 0, unique_conic_forms)
    cache_conic_form!(unique_conic_forms, x, objective)
  end
  return get_conic_form(unique_conic_forms, x)
end

sum_largest(x::AbstractExpr, k::Int) = SumLargestAtom(x, k)
sum_smallest(x::AbstractExpr, k::Int) = -SumLargestAtom(-x, k)

# This atom computes sort(w)'*sort(x)
# ie if w and x are sorted in decreasing order, then we return w'*x
# If w = [1 1 1 0 0 0 ... 0], it computes the sum of the 3 largest elements of x
type WeightedSumAtom <: AbstractExpr
  head::Symbol
  id_hash::Uint64
  children::(AbstractExpr,)
  size::(Int, Int)
  w::Vector

  function WeightedSumAtom(x::AbstractExpr, w::Vector)
    if !(length(w) == get_vectorized_size(x))
      error("x and w must be the same size")
    end
    children = (x,)
    return new(:weighted_sum, hash((children, w)), children, (1,1), w)
  end
end

function sign(x::WeightedSumAtom)
  if all(x.w.>=0)
    return sign(x.children[1])
  elseif all(x.w.<=0)
    return sign(x.children[1])
  else
    return NoSign()
  end
end

function monotonicity(x::WeightedSumAtom)
  if all(x.w.>=0)
    return (Nondecreasing(), )
  else
    return (NoMonotonicity(), )
  end
end

function curvature(x::WeightedSumAtom)
  return ConvexVexity()
end

function evaluate(x::WeightedSumAtom)
  return sum(sort(vec(evaluate(x.children[1])), rev=true) .* sort(vec(x.w), rev=true))
end

function conic_form!(x::WeightedSumAtom, unique_conic_forms::UniqueConicForms)
  if !has_conic_form(unique_conic_forms, x)
    y = x.children[1]
    if all(size(y)) > 1
      y = vec(y)
      w = vec(x.w)
    else
      w = x.w
    end
    mu = Variable(size(y))
    nu = Variable(size(y))
    onesvec = ones(size(y))
    # given by the solution to
    # minimize sum(mu) + sum(nu)
    # subject to y*w' <= onesvec*nu' + mu*onesvec'
    objective = conic_form!(sum(mu) + sum(nu), unique_conic_forms)
    conic_form!(y*w' <= onesvec*nu' + mu*onesvec', unique_conic_forms)
    cache_conic_form!(unique_conic_forms, x, objective)
  end
  return get_conic_form(unique_conic_forms, x)
end
weighted_sum(x::AbstractExpr, w::Vector) = WeightedSumAtom(x, w)
