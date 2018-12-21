#############################################################################
# conv.jl
# Handles convolution between a constant vector and an expression vector.
#############################################################################

export conv

function conv(x::Value, y::AbstractExpr)
    if (size(x,2) != 1 && length(size(x)) != 1) || size(y,2) != 1
        error("convolution only supported between two vectors")
    end
    m = length(x)
    n = y.size[1]
    X = zeros(m+n - 1, n)
    for i = 1:n
        X[i:m+i-1, i] = x
    end
    return X*y
end

conv(x::AbstractExpr, y::Value) = conv(y, x)
