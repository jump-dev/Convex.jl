export inner_product

function inner_product(x::AbstractExpr, y::AbstractExpr)
    if x.size == y.size && x.size[1] == x.size[2]
        return real(tr(x'*y))
    else
        error("arguments must be square matrix of same dimension")
    end
end

inner_product(x::Value, y::AbstractExpr) = inner_product(Constant(x),y)
inner_product(x::AbstractExpr, y::Value) = inner_product(x,Constant(y))
