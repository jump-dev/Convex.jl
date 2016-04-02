export geomean_eig, nthroot_det

## to do: need a geomean function that takes vectors of length n
# and computes the nth root of the product of their entries
# we only have this for n=2

type NthRootDetAtom <: AbstractExpr
  head::Symbol
  id_hash::UInt64
  children::Tuple{AbstractExpr}
  size::Tuple{Int, Int}

  function NthRootDetAtom(x::AbstractExpr)
    children = (x,)
    return new(:nthrootdet, hash(children), children, (1, 1))
  end
end

function sign(x::NthRootDetAtom)
  return Pos()
end

function monotonicity(x::NthRootDetAtom)
  return (NoMonotonicity(),)
end

function curvature(x::NthRootDetAtom)
  return ConcaveVexity()
end

function evaluate(x::NthRootDetAtom)
  n = size(x,1)
  return (det(evaluate(x.children[1])))^(1/n)
end

function conic_form!(x::NthRootDetAtom, unique_conic_forms::UniqueConicForms)
  if !has_conic_form(unique_conic_forms, x)
    A = x.children[1]
    D = Variable(size(A)) # diagonal matrix
    U = Variable(size(A)) # upper triangular matrix

    # objective given by the geometric mean of the eigenvalues,
    # represented by the diagonal of diagonal matrix D
    objective = conic_form!(geomean(diag(D)), unique_conic_forms)

    # force D to be diagonal; for U to be upper triangular
    for i in 1:A.size[1]
      for j in 1:A.size[2]
        if i != j
          conic_form!(D[i,j] == 0, unique_conic_forms)
        end
        if i > j
          conic_form!(U[i,j] == 0, unique_conic_forms)
        end
      end
    end

    # diagonals of D and U are the same
    conic_form!(diag(D) == diag(U), unique_conic_forms)

    # A and [D U; U' A] need to be positive semidefinite
    conic_form!(isposdef(A), unique_conic_forms)
    conic_form!(isposdef([D U; U' A]), unique_conic_forms)

    cache_conic_form!(unique_conic_forms, x, objective)
  end
  return get_conic_form(unique_conic_forms, x)
end

geomean_eig(x::AbstractExpr) = NthRootDetAtom(x)
nthroot_det(x::AbstractExpr) = NthRootDetAtom(x)
