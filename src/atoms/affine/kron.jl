import Base.kron
export kron


type KronAtom <: AbstractExpr
  head::Symbol
  id_hash::UInt64
  children::Tuple{AbstractExpr, AbstractExpr}
  size::Tuple{Int, Int}

  function KronAtom(x::AbstractExpr, y::AbstractExpr)
    if vexity(x.children[1]) != ConstVexity() && vexity(x.children[1]) != ConstVexity()
        error("Kron of two non-constant expressions is not DCP compliant")
    else
        sz = (size(x)[1]*size(y)[1], size(x)[2]*size(y)[2])
        children = (x, y)
        return new(:kron, hash(children), children, sz)
    end
end

function sign(x::KronAtom)
  return sign(x.children[1]) * sign(x.children[2])
end

function monotonicity(x::KronAtom)
  return (sign(x.children[2]) * Nondecreasing(), sign(x.children[1]) * Nondecreasing())
end


function curvature(x::KronAtom)
    return ConstVexity()
  end
end

function evaluate(x::KronAtom)
  return kron(evaluate(x.children[1]),evaluate(x.children[2]))
end

function conic_form!(x::KronAtom, unique_conic_forms::UniqueConicForms)
  if !has_conic_form(unique_conic_forms, x)
    # scalar multiplication
    if x.children[1].size == (1, 1) || x.children[2].size == (1, 1)
      if vexity(x.children[1]) == ConstVexity()
        const_child = x.children[1]
        expr_child = x.children[2]
      elseif vexity(x.children[2]) == ConstVexity()
        const_child = x.children[2]
        expr_child = x.children[1]
      else
        error("multiplication of two non-constant expressions is not DCP compliant")
      end
      objective = conic_form!(expr_child, unique_conic_forms)

      # make sure all 1x1 sized objects are interpreted as scalars, since
      # [1] * [1, 2, 3] is illegal in julia, but 1 * [1, 2, 3] is ok
      if const_child.size == (1, 1)
        const_multiplier = evaluate(const_child)[1]
      else
        const_multiplier = reshape(evaluate(const_child), get_vectorized_size(const_child), 1)
      end

      objective = const_multiplier * objective

    # left matrix multiplication
    elseif x.children[1].head == :constant
      objective = conic_form!(x.children[2], unique_conic_forms)
      objective = kron(speye(x.size[2]), x.children[1].value) * objective
    # right matrix multiplication
    else
      objective = conic_form!(x.children[1], unique_conic_forms)
      objective = kron(x.children[2].value', speye(x.size[1])) * objective
    end
    cache_conic_form!(unique_conic_forms, x, objective)
  end
  return get_conic_form(unique_conic_forms, x)
end



