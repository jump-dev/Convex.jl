import Base.kron
export kron
export sign, monotonicity, curvature, evaluate, conic_form!

type KronAtom <: AbstractExpr
  head::Symbol
  id_hash::UInt64
  children::Tuple{AbstractExpr, AbstractExpr}
  size::Tuple{Int, Int}

  function KronAtom(x::AbstractExpr, y::AbstractExpr)
    if vexity(x) != ConstVexity() && vexity(y) != ConstVexity()
      error("Kron of two non-constant expressions is not DCP compliant")
    else
      sz = (size(x)[1]*size(y)[1], size(x)[2]*size(y)[2])
      children = (x, y)
      return new(:kron, hash(children), children, sz)
    end
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

function evaluate(x::KronAtom)
  return kron(evaluate(x.children[1]),evaluate(x.children[2]))
end

function conic_form!(x::KronAtom, unique_conic_forms::UniqueConicForms)
  if !has_conic_form(unique_conic_forms, x)
    # scalar multiplication
    # if x.children[1].size == (1, 1) || x.children[2].size == (1, 1)
    #   if vexity(x.children[1]) == ConstVexity()
    #     const_child = x.children[1]
    #     expr_child = x.children[2]
    #   elseif vexity(x.children[2]) == ConstVexity()
    #     const_child = x.children[2]
    #     expr_child = x.children[1]
    #   else
    #     error("multiplication of two non-constant expressions is not DCP compliant")
    #   end
    #   objective = conic_form!(expr_child, unique_conic_forms)

    #   # make sure all 1x1 sized objects are interpreted as scalars, since
    #   # [1] * [1, 2, 3] is illegal in julia, but 1 * [1, 2, 3] is ok
    #   if const_child.size == (1, 1)
    #     const_multiplier = evaluate(const_child)[1]
    #   else
    #     const_multiplier = reshape(evaluate(const_child), get_vectorized_size(const_child), 1)
    #   end

    #   objective = const_multiplier * objective

    # # left matrix multiplication
    # else x.children[1].head == :constant
      objective = conic_form!(x.children[2], unique_conic_forms)
      a = evaluate(x.children[1])
      for key in objective.keys
        rows1 = SparseMatrixCSC{Float64,Int32}
        rows2 = SparseMatrixCSC{Float64,Int32}
        for i in 1:size(x.children[1])[1]
          row1 = SparseMatrixCSC{Float64,Int32}
          rows2 = SparseMatrixCSC{Float64,Int32}
          for j in 1:size(x.children[1])[2]
            x = objective[key][1]*a[i,j]
            y = objective[key][2]*a[i,j]
            push!(row1,x)
            push!(row2,y)
          end
          push!(rows1, foldl(hcat, row1))
          push!(rows2, foldl(hcat, row2))
        end
        objective[key][1] = foldl(vcat, rows1)
        objective[key][2] = foldl(vcat, rows2)
      end

      #objective = kron(x.children[1].value,ones(x.size[2])) * objective
    # right matrix multiplication
    #else   when second argumet is constant
     # objective = conic_form!(x.children[1], unique_conic_forms)
      #objective = kron(x.children[2].value', speye(x.size[1])) * objective
    #end
    cache_conic_form!(unique_conic_forms, x, objective)
  end
  return get_conic_form(unique_conic_forms, x)
end

kron(x::Value, y::AbstractExpr) = KronAtom(Constant(x), y)

