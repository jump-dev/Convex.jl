function conv(x::Value, y::AbstractExpr)
    if length(x) != size(x, 1) || size(y, 2) > 1
        error("convolution only supported between two vectors")
    end
    m, n = length(x), size(y, 1)
    X = spzeros(eltype(x), m + n - 1, n)
    for i in 1:n, j in 1:m
        X[i+j-1, i] = x[j]
    end
    return X * y
end

conv(x::AbstractExpr, y::Value) = conv(y, x)
