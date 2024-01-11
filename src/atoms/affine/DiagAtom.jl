mutable struct DiagAtom <: AbstractExpr
    children::Tuple{AbstractExpr}
    size::Tuple{Int,Int}
    k::Int

    function DiagAtom(x::AbstractExpr, k::Int = 0)
        K = min(x.size[1], x.size[2])
        if !(-K < k < K)
            error("Bounds error in calling diag")
        end
        return new((x,), (K - k, 1), k)
    end
end

head(io::IO, ::DiagAtom) = print(io, "diag")

Base.sign(x::DiagAtom) = sign(x.children[1])

monotonicity(::DiagAtom) = (Nondecreasing(),)

curvature(::DiagAtom) = ConstVexity()

function evaluate(x::DiagAtom)
    return LinearAlgebra.diag(evaluate(x.children[1]), x.k)
end

LinearAlgebra.diag(x::AbstractExpr, k::Int = 0) = DiagAtom(x, k)

# Finds the "k"-th diagonal of x as a column vector.
#
# If k == 0, it returns the main diagonal and so on.
#
# Let x be of size m x n and d be the diagonal.
#
# Since x is vectorized, the way canonicalization works is:
#
# 1. We calculate the size of the diagonal (sz_diag) and the first index
#    of vectorized x that will be part of d
# 2. We create the coefficient matrix for vectorized x, called coeff of size
#    sz_diag x mn
# 3. We populate coeff with 1s at the correct indices
#
# The canonical form will then be: coeff * x - d = 0
function new_conic_form!(context::Context{T}, x::DiagAtom) where {T}
    num_rows, num_cols = x.children[1].size
    if x.k >= 0
        start_index = x.k * num_rows + 1
        sz_diag = Base.min(num_rows, num_cols - x.k)
    else
        start_index = -x.k + 1
        sz_diag = Base.min(num_rows + x.k, num_cols)
    end
    select_diag = spzeros(T, sz_diag, length(x.children[1]))
    for i in 1:sz_diag
        select_diag[i, start_index] = 1
        start_index += num_rows + 1
    end
    child_obj = conic_form!(context, only(AbstractTrees.children(x)))
    return operate(add_operation, T, sign(x), select_diag, child_obj)
end
