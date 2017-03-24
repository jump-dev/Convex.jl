import Base.kron
export kron


#### TODO: wite the conic_form implemenatatioo

function kron(a::Value, b::Convex.Variable)
  rows = Convex.AbstractExpr[]
  for i in 1:size(a)[1]
    row = Convex.AbstractExpr[]
    for j in 1:size(a)[2]
      if isreal(a[i,j])
        push!(row, real(a[i, j]) * b)
      else
        push!(row, a[i, j] * b)
      end        
    end
    push!(rows, foldl(hcat, row))
  end
  return foldl(vcat, rows)
end


function kron(a::Convex.Variable, b::Union{AbstractArray, Convex.Constant})
  rows = Convex.AbstractExpr[]
  b = Constant(b)
  for i in 1:size(a)[1]
    row = Convex.AbstractExpr[]
    for j in 1:size(a)[2]
      push!(row, a[i, j] * b)
    end
    push!(rows, foldl(hcat, row))
  end
  return foldl(vcat, rows)
end
