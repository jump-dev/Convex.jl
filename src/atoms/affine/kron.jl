function LinearAlgebra.kron(a::Value, b::AbstractExpr)
    rows = AbstractExpr[]
    for i in 1:size(a)[1]
        row = AbstractExpr[]
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


function LinearAlgebra.kron(a::AbstractExpr, b::Value)
    rows = AbstractExpr[]
    for i in 1:size(a)[1]
        row = AbstractExpr[]
        for j in 1:size(a)[2]
            if isreal(b)
                push!(row, a[i, j] * real(b))
            else
                push!(row, a[i, j] * b)
            end
        end
        push!(rows, foldl(hcat, row))
    end
    return foldl(vcat, rows)
end
