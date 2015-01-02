export logdet

type LogDetAtom <: AbstractExpr
  head::Symbol
  id_hash::Uint64
  children::(AbstractExpr,)
  size::(Int, Int)

  function LogDetAtom(x::AbstractExpr)
    children = (x,)
    return new(:logdet, hash(children), children, (1, 1))
  end
end

function sign(x::LogDetAtom)
  return NoSign()
end

function monotonicity(x::LogDetAtom)
  return (NoMonotonicity(),)
end

function curvature(x::LogDetAtom)
  return ConcaveVexity()
end

function evaluate(x::LogDetAtom)
  return log(det(evaluate(x.children[1])))
end

function conic_form!(x::LogDetAtom, unique_conic_forms::UniqueConicForms)
  if !has_conic_form(unique_conic_forms, x)
    A = x.children[1]
    D = Variable(size(A)) # diagonal matrix
    U = Variable(size(A)) # upper triangular matrix

    # objective given by the sum of the log of diagonal matrix D
    objective = conic_form!(sum(log(diag(D))), unique_conic_forms)

    # enforce D and U to be diagonal and upper triangular
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

logdet(x::AbstractExpr) = LogDetAtom(x)
