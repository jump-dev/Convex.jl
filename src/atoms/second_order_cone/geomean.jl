import Base.sqrt
export GeoMeanAtom, geomean, sqrt, sgeomean
export sign, monotonicity, curvature, conic_form!

type GeoMeanAtom <: AbstractExpr
  head::Symbol
  id_hash::UInt64
  children::Tuple{AbstractExpr, AbstractExpr}
  size::Tuple{Int, Int}

  function GeoMeanAtom(x::AbstractExpr, y::AbstractExpr)
    if x.size != y.size
      error("geo mean must take two arguments of the same size")
    end
    children = (x, y)
    return new(:geomean, hash(children), children, x.size)
  end
end

function sign(q::GeoMeanAtom)
  return Positive()
end

function monotonicity(q::GeoMeanAtom)
  return (Nondecreasing(), Nondecreasing())
end

function curvature(q::GeoMeanAtom)
  return ConcaveVexity()
end

function evaluate(q::GeoMeanAtom)
  return sqrt(evaluate(q.children[1]) .* evaluate(q.children[2]))
end

function conic_form!(q::GeoMeanAtom, unique_conic_forms::UniqueConicForms)
  if !has_conic_form(unique_conic_forms, q)
    sz = q.children[1].size
    t = Variable(sz[1], sz[2])
    qol_objective = conic_form!(t, unique_conic_forms)
    x, y = q.children
    conic_form!(SOCElemConstraint(y + x, y - x, 2 * t), unique_conic_forms)
    conic_form!(x >= 0, unique_conic_forms)
    conic_form!(y >= 0, unique_conic_forms)
    cache_conic_form!(unique_conic_forms, q, qol_objective)
  end
  return get_conic_form(unique_conic_forms, q)
end

geomean(x::AbstractExpr, y::AbstractExpr) = GeoMeanAtom(x, y)
sqrt(x::AbstractExpr) = GeoMeanAtom(x, Constant(ones(x.size[1], x.size[2])))

# (almost) extend geometric mean to vectors of length n via recursion
# geomean(x) = prod(x)^(1/\bar n)
# where \bar n is the smallest power of 2 bigger than n

function power_of_2_gt(n::Int)
  Int(2^ceil(log(2,n)))
end

function sgeomean(x::AbstractExpr)
  if length(x) > 2
    nbar = power_of_2_gt(length(x))
    half_nbar = Int(nbar/2)
    first_half = x[1:half_nbar]
    # append ones to last half until it's a power of 2
    last_half = vcat(vec(x[half_nbar+1:end]), ones(nbar - length(x)))
    return geomean(sgeomean(first_half), 
                 sgeomean(last_half))
  elseif length(x) == 2  
    return geomean(x[1],x[2])
  else
    return x
  end
end

# a test for sgeomean
# p = maximize(sgeomean(x), x<= 2, x[1]<=1); solve!(p);
# @assert p.optval == prod(x.value)^(1/16)