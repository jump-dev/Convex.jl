import Base.reshape, Base.vec
export reshape, vec

# Creates a new expression, same as x, but updates the size
function reshape(x::AbstractCvxExpr, m::Int64, n::Int64)
  sz = get_vectorized_size(x.size)
  if m * n != sz
    error("Can't reshape expression with $sz entries into $(m * n)")
  end

  this = CvxExpr(:reshape, [x], x.vexity, x.sign, (m, n))
  coeffs = VecOrMatOrSparse[speye(sz), -speye(sz)]
  vars = [x.uid, this.uid]
  canon_constr = CanonicalConstr(coeffs, vars, spzeros(sz, 1), true, false)

  canon_constr_array = x.canon_form()
  push!(canon_constr_array, canon_constr)

  this.canon_form = ()->canon_constr_array
  this.evaluate = ()->Base.reshape(x.evaluate(), m, n)
  return this
end

function vec(x::AbstractCvxExpr)
  return reshape(x, get_vectorized_size(x), 1)
end
