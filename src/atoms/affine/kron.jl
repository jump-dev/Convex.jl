import Base.kron
export kron

function kron(a::Union{AbstractArray, Convex.Constant}, b::Convex.Variable)
  rows = Convex.AbstractExpr[]
  for i in 1:size(a)[1]
    row = Convex.AbstractExpr[]
    for j in 1:size(a)[2]
      push!(row, a[i, j] * b)
    end
    push!(rows, foldl(hcat, row))
  end
  return foldl(vcat, rows)
end


function kron(a::Convex.Variable, b::Union{AbstractArray, Convex.Constant})
  rows = Convex.AbstractExpr[]
  for i in 1:size(b)[1]
    row = Convex.AbstractExpr[]
    for j in 1:size(b)[2]
      push!(row, b[i, j] * a)
    end
    push!(rows, foldl(vcat, row))
  end
  return foldl(hcat, rows)
end