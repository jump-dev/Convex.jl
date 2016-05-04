import Base: eig
export eig

# schatten(A) or eig(A) computes the eigenvalue of a semidefinite matrix A

type SchattenAtom <: AbstractExpr
  head::Symbol
  id_hash::UInt64
  children::Tuple{AbstractExpr}
  size::Tuple{Int, Int}

  function SchattenAtom(x::AbstractExpr)
    children = (x,)
    if size(x,1) != size(x,2)
      error("schatten can only be applied to square matrices")
    end
    return new(:schatten, hash(children), children, (size(x,1), 1))
  end
end

function sign(x::SchattenAtom)
  return Pos()
end

function monotonicity(x::SchattenAtom)
  return (NoMonotonicity(),)
end

function curvature(x::SchattenAtom)
  return AffineVexity() # XXX does this make any sense?
end

function evaluate(x::SchattenAtom)
  A = evaluate(x.children[1])
  l,v = eig(A)
  return l
end

function conic_form!(x::SchattenAtom, unique_conic_forms::UniqueConicForms)
  if !has_conic_form(unique_conic_forms, x)
    A = x.children[1]
    d = Variable(size(A,1)) # vector
    U = Variable(size(A)) # upper triangular matrix

    # objective given by the function f of the eigenvalues,
    # represented by the diagonal of diagonal matrix D
    objective = conic_form!(d, unique_conic_forms)

    # force U to be upper triangular
    for i in 1:A.size[1]
      for j in 1:A.size[2]
        if i > j
          conic_form!(U[i,j] == 0, unique_conic_forms)
        end
      end
    end

    # diagonals of D and U are the same
    conic_form!(d == diag(U), unique_conic_forms)

    # A and [D U; U' A] need to be positive semidefinite
    conic_form!(isposdef(A), unique_conic_forms)
    conic_form!(isposdef([diagm(d) U; U' A]), unique_conic_forms)

    cache_conic_form!(unique_conic_forms, x, objective)
  end
  return get_conic_form(unique_conic_forms, x)
end

eig(x::AbstractExpr) = SchattenAtom(x)
