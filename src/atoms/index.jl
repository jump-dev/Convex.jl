export getindex

# TODO: Only works for vectors
function getindex(x::AbstractCvxExpr, indices::Range1{Int64})
  c = spzeros(x.size...)
  #TODO: Does not work for range @kvmohan
  c[indices, :] = 1.0
  return dot(c, x)
end

getindex(x::AbstractCvxExpr, index::Int64) = getindex(x, index:index)
# If a non-Int64 is passed, we try to convert it into an Int64. For something like 1.0
# this will work. For a float like 1.5, convert will throw an InexactError
getindex(x::AbstractCvxExpr, index::Number) = getindex(x, convert(Int64, index))
