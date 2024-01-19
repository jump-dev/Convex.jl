# How to create your own atom
# This file is a good example of how a user can create his or her own atom
# This antidiag atom is then included in n_queens.jl to solve the n queens problem

#############################################################################
# antidiag.jl
# Returns the kth anti-diagonal (counterdiagonal, secondary diagonal, or minor diagonal) of a matrix expression
# All expressions and atoms are subtpyes of AbstractExpr.
# Please read expressions.jl first.
#############################################################################

import Convex
export antidiag

### Diagonal
### Represents the kth diagonal of an mxn matrix as a (min(m, n) - k) x 1 vector
struct AntidiagAtom <: Convex.AbstractExpr
    head::Symbol
    id_hash::UInt64
    children::Tuple{Convex.AbstractExpr}
    size::Tuple{Int,Int}
    k::Int

    function AntidiagAtom(x::Convex.AbstractExpr, k::Int = 0)
        num_rows, num_cols = x.size
        if !(-num_rows <= k <= num_cols)
            error("Bounds error in calling diag")
        end
        children = (x,)
        return new(
            :antidiag,
            hash((children, k)),
            children,
            (minimum(x.size) - k, 1),
            k,
        )
    end
end

Base.sign(x::AntidiagAtom) = sign(x.children[1])

Convex.monotonicity(::AntidiagAtom) = (Convex.Nondecreasing(),)

Convex.curvature(::AntidiagAtom) = Convex.ConstVexity()

function Convex.evaluate(x::AntidiagAtom)
    return diag(reverse(Convex.evaluate(x.children[1]), dims = 1), x.k)
end

antidiag(x::Convex.AbstractExpr, k::Int = 0) = AntidiagAtom(x, k)

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
function Convex.new_conic_form!(
    context::Convex.Context{T},
    x::AntidiagAtom,
) where {T}
    (num_rows, num_cols) = x.children[1].size
    k = x.k
    if k >= 0
        start_index = k * num_rows + num_rows
        sz_diag = Base.min(num_rows, num_cols - k)
    else
        start_index = num_rows + k
        sz_diag = Base.min(num_rows + k, num_cols)
    end
    select_diag = spzeros(T, sz_diag, length(x.children[1]))
    for i in 1:sz_diag
        select_diag[i, start_index] = 1
        start_index += num_rows - 1
    end
    objective = Convex.conic_form!(context, only(Convex.children(x)))
    return Convex.operate(
        Convex.add_operation,
        T,
        Convex.sign(x),
        select_diag,
        objective,
    )
end
