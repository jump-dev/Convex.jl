export getindex

# TODO:
# Allow indexing with just `:`. eg, x[1, :] or x[2, 5:]

# Slice a variable or an expression to get only certain elements
# The canonical form is straightforward:
# 1. Calculate number of elements needed, call this num_e
# 2. Create coefficient matrix, coeff, of size num_e * get_vectorized_size(x)
# 3. Add 1's at the correct places of coeff to allow
# coeff * vectorized(x) - vectorized(indexed_elements) = 0
function getindex(x::AbstractCvxExpr, rows::AbstractArray, cols::AbstractArray=[1])
  if size(rows, 1) != 1 && size(rows, 2) != 1
    error("Expected a vector but got size $(size(rows))")
  end
  if size(cols, 1) != 1 && size(cols, 2) != 1
    error("Expected a vector but got size $(size(cols))")
  end

  # Get rid of duplicate elements
  rows = unique(rows)
  cols = unique(cols)

  # Number of rows/cols in the coefficient for x in our canonical form
  num_rows_coeff = length(rows) * length(cols)
  num_cols_coeff = get_vectorized_size(x)

  coeff = spzeros(num_rows_coeff, num_cols_coeff)

  # Create the coeff matrix such that coeff * vec(x) = vec(x[rows, cols])
  k = 1
  num_rows = x.size[1]
  for c in cols
    for r in rows
      idx = num_rows * (c - 1) + r
      coeff[k, idx] = 1
      k += 1
    end
  end

  this = CvxExpr(:index, [x], x.vexity, x.sign, (length(rows), length(cols)))
  coeffs = VecOrMatOrSparse[coeff, -speye(num_rows_coeff)]
  vars = [x.uid, this.uid]
  constant = zeros(num_rows_coeff)

  canon_constr_array = [CanonicalConstr(coeffs, vars, constant, true, false)]
  append!(canon_constr_array, x.canon_form())
  this.canon_form = ()->canon_constr_array

  this.evaluate = ()->x.evaluate()[rows, cols]

  return this
end

getindex(x::AbstractCvxExpr, row::Int64, col::Int64=1) = getindex(x, row:row, col:col)
getindex(x::AbstractCvxExpr, row::Number, col::Number=1) = getindex(x, row:row, col:col)
getindex(x::AbstractCvxExpr, row::Int64, cols::AbstractArray) = getindex(x, row:row, cols)
getindex(x::AbstractCvxExpr, rows::AbstractArray, col::Int64) = getindex(x, rows, col:col)

# If a non-Int64 is passed, we try to convert it into an Int64. For something like 1.0
# this will work. For a float like 1.5, convert will throw an InexactError
# Note that we don't use int(index) because `int` rounds the number whereas float
# wil only convert floats of the form 5.0, 4.0 etc. and not 4.5
getindex(x::AbstractCvxExpr, index::Number) = getindex(x, convert(Int64, index))
