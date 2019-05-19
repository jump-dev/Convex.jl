#############################################################################
# diag.jl
# Returns the kth diagonal of a matrix expression
# All expressions and atoms are subtpyes of AbstractExpr.
# Please read expressions.jl first.
#############################################################################

# k >= min(num_cols, num_rows) || k <= -min(num_rows, num_cols)
import LinearAlgebra.diag
export diag
#export sign, curvature, monotonicity, evaluate

### Diagonal
### Represents the kth diagonal of an mxn matrix as a (min(m, n) - k) x 1 vector

struct DiagAtom <: AbstractExpr
    head::Symbol
    id_hash::UInt64
    children::Tuple{AbstractExpr}
    size::Tuple{Int, Int}
    k::Int

    function DiagAtom(x::AbstractExpr, k::Int=0)
        (num_rows, num_cols) = x.size

        if k >= min(num_rows, num_cols) || k <= -min(num_rows, num_cols)
            error("bounds error in calling diag")
        end

        children = (x, )
        return new(:diag, hash((children, k)), children, (min(num_rows, num_cols) - k, 1), k)
    end
end

## Type Definition Ends


function sign(x::DiagAtom)
    return sign(x.children[1])
end

# The monotonicity
function monotonicity(x::DiagAtom)
    return (Nondecreasing(),)
end

# If we have h(x) = f o g(x), the chain rule says h''(x) = g'(x)^T f''(g(x))g'(x) + f'(g(x))g''(x);
# this represents the first term
function curvature(x::DiagAtom)
    return ConstVexity()
end

function evaluate(x::DiagAtom)
    return diag(evaluate(x.children[1]), x.k)
end

## API begins
diag(x::AbstractExpr, k::Int=0) = DiagAtom(x, k)
## API ends

# Finds the "k"-th diagonal of x as a column vector
# If k == 0, it returns the main diagonal and so on
# Let x be of size m x n and d be the diagonal
# Since x is vectorized, the way canonicalization works is:
#
# 1. We calculate the size of the diagonal (sz_diag) and the first index
# of vectorized x that will be part of d
# 2. We create the coefficient matrix for vectorized x, called coeff of size
# sz_diag x mn
# 3. We populate coeff with 1s at the correct indices
# The canonical form will then be:
# coeff * x - d = 0
function conic_form!(x::DiagAtom, unique_conic_forms::UniqueConicForms=UniqueConicForms())
    if !has_conic_form(unique_conic_forms, x)
        (num_rows, num_cols) = x.children[1].size
        k = x.k

        if k >= 0
            start_index = k * num_rows + 1
            sz_diag = Base.min(num_rows, num_cols - k)
        else
            start_index = -k + 1
            sz_diag = Base.min(num_rows + k, num_cols)
        end

        select_diag = spzeros(sz_diag, length(x.children[1]))
        for i in 1:sz_diag
            select_diag[i, start_index] = 1
            start_index += num_rows + 1
        end

        objective = conic_form!(x.children[1], unique_conic_forms)
        new_obj = select_diag * objective
        cache_conic_form!(unique_conic_forms, x, new_obj)
    end
    return get_conic_form(unique_conic_forms, x)
end
