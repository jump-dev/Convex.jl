export getindex

function getindex(x::AbstractCvxExpr, rows::Range1{Int64}, cols::Range1{Int64}=1:1)
  # Overload for row vectors
  if x.size[1] == 1 && rows != 1:1
    cols = rows
    rows = 1:1
  end

  length_rows = length(rows)
  lh = spzeros(length_rows, x.size[1])
  lh[:, rows] = speye(length_rows)

  length_cols = length(cols)
  rh = spzeros(x.size[2], length_cols)
  rh[cols, :] = speye(length_cols)

  return lh * x * rh
end

getindex(x::AbstractCvxExpr, row::Int64, col::Int64=1) = getindex(x, row:row, col:col)

# If a non-Int64 is passed, we try to convert it into an Int64. For something like 1.0
# this will work. For a float like 1.5, convert will throw an InexactError
getindex(x::AbstractCvxExpr, index::Number) = getindex(x, convert(Int64, index))
