function quadform(x::Value, A::AbstractExpr)
    return x' * A * x
end

is_psd(A::AbstractMatrix{T}) where {T} = isposdef(A + eps(float(T))*Matrix{T}(I, size(A)...))

function quadform(x::AbstractExpr, A::Value)
    if length(size(A)) != 2 || size(A, 1) != size(A, 2)
        error("Quadratic form only takes square matrices")
    end
    if !ishermitian(A)
        error("Quadratic form only defined for Hermitian matrices")
    end
    if is_psd(A)
        factor = 1
    elseif is_psd(-A)
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
