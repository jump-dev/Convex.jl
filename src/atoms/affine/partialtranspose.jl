import Base.sign
export partialtranspose

function partialtranspose(x::AbstractExpr, sys::Int, dim::Vector)
    if x.size[1] != x.size[2]
        error("Only square matrices are supported")
    end
    if ! (1 ≤ sys ≤ length(dims))
        error("Invalid system, should between 1 and ", length(dims), " got ", sys)
    end
    if x.size[1] ≠ prod(dims)
        error("Dimension of system doesn't correspond to dimension of subsystems")
    end
    if !is_square(x.size[1])
        error("Size of the matrix must be a perfect square")
    end

    



end

function is_square(Int::a)
    root = sqrt(a)
    if root*root == a
        return true
    else
        return false
    end
end



