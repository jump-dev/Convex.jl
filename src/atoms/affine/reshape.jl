import Base.reshape, Base.vec
export reshape, vec
export sign, curvature, monotonicity, evaluate, conic_form


type ReshapeAtom <: AbstractExpr
  head::Symbol
  id::Uint64
  children::(AbstractExpr,)
  size::(Int64, Int64)

  function ReshapeAtom(x::AbstractExpr, m::Int64, n::Int64)
    if m * n != get_vectorized_size(x)
      error("Cannot reshape expression of size $(x.size) to ($(m), $(n))")
    end
    return new(:reshape, object_id(x), (x,), (m, n))
  end
end

function sign(x::ReshapeAtom)
  return sign(x.children[1])
end

function monotonicity(x::ReshapeAtom)
  return (Nondecreasing(),)
end

function curvature(x::ReshapeAtom)
  return ConstVexity()
end

function evaluate(x::ReshapeAtom)
  return reshape(evaluate(x.children[1]), x.size[1], x.size[2])
end

function conic_form(x::ReshapeAtom, unique_constr)
  conic_form(x.children[1], unique_constr)
end


reshape(x::AbstractExpr, m::Int64, n::Int64) = ReshapeAtom(x, m, n)
vec(x::AbstractExpr) = reshape(x, get_vectorized_size(x), 1)
