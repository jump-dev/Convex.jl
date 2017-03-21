import Base.kron
export kron

function kron(a::Union{Value, Constant}, b::AbstractExpr)
  rows = AbstractExpr[]
  a = Constant(a)
  for i in 1:size(a)[1]
    row = AbstractExpr[]
    for j in 1:size(a)[2]
      push!(row, a[i, j] * b)
    end
    push!(rows, foldl(hcat, row))
  end
  return foldl(vcat, rows)
end


function kron(a::AbstractExpr, b::Union{Value, Constant})
  rows = AbstractExpr[]
  b = Constant(b)
  for i in 1:size(a)[1]
    row = AbstractExpr[]
    for j in 1:size(a)[2]
      push!(row, a[i, j] * b)
    end
    push!(rows, foldl(hcat, row))
  end
  return foldl(vcat, rows)
end

