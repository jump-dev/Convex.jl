function inner_product(x::AbstractExpr, y::AbstractExpr)
    if x.size == y.size && x.size[1] == x.size[2]
        return real(tr(x'*y))
    else
        error("Arguments must be square matrix of same dimension")
    end
end

inner_product(x::Value, y::AbstractExpr) = inner_product(constant(x),y)
inner_product(x::AbstractExpr, y::Value) = inner_product(x,constant(y))
