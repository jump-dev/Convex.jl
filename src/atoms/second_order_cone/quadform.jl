function quadform(x::Value, A::AbstractExpr)
    return x' * A * x
end

function quadform(x::AbstractExpr, A::Value)
    if length(size(A)) != 2 || size(A, 1) != size(A, 2)
        error("Quadratic form only takes square matrices")
    end
    if !ishermitian(A)
        error("Quadratic form only defined for Hermitian matrices")
    end
    if isposdef(A)
        factor = 1
    elseif isposdef(.-A)
        factor = -1
    else
        error("Quadratic forms supported only for semidefinite matrices")
    end

    P = sqrt(Hermitian(factor * A))
    return factor * square(norm2(P * x))
end

function quadform(x::AbstractExpr, A::AbstractExpr)
    if vexity(x) == ConstVexity()
        return quadform(evaluate(x), A)
    elseif vexity(A) == ConstVexity()
        return quadform(x, evaluate(A))
    else
        error("Either `x` or `A` must be constant in `quadform(x,A)`.")
    end
end
