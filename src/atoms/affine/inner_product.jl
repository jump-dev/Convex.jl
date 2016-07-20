export inner_product

function inner_product(x::AbstractExpr,y::AbstractExpr)
    if x.size==y.size
        if x.size[1] == x.size[2]
            if y.size[1] == y.size[2]
                return(trace(y*x))
            else 
                error("Second argument must be a square matrix")
        else 
            error("First argument must be a square matrix")

    else 
        error("Dimension of both the matrix must be same")
end

#inner_product(x::AbstractExpr, y::AbstractExpr) = transpose(x) * y
inner_product(x::Value, y::AbstractExpr) = inner_product(Constant(x),y)
inner_product(x::AbstractExpr, y::Value) = inner_product(x,Constant(y))