#############################################################################
# diag.jl
# Returns the kth diagonal of a matrix expression
# All expressions and atoms are subtpyes of AbstractExpr.
# Please read expressions.jl first.
#############################################################################

# k >= min(num_cols, num_rows) || k <= -min(num_rows, num_cols)

### Diagonal
### Represents the kth diagonal of an mxn matrix as a (min(m, n) - k) x 1 vector

struct DiagAtom <: AbstractExpr
    children::Tuple{AbstractExpr}
    size::Tuple{Int,Int}
    k::Int

    function DiagAtom(x::AbstractExpr, k::Int = 0)
        (num_rows, num_cols) = x.size

        if k >= min(num_rows, num_cols) || k <= -min(num_rows, num_cols)
            error("Bounds error in calling diag")
        end

        children = (x,)
        return new(children, (min(num_rows, num_cols) - k, 1), k)
    end
end

head(io::IO, ::DiagAtom) = print(io, "diag")

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
LinearAlgebra.diag(x::AbstractExpr, k::Int = 0) = DiagAtom(x, k)
## API ends

function conic_form!(context::Context{T}, x::DiagAtom) where {T}
    (num_rows, num_cols) = x.children[1].size
    k = x.k

    if k >= 0
        start_index = k * num_rows + 1
        sz_diag = Base.min(num_rows, num_cols - k)
    else
        start_index = -k + 1
        sz_diag = Base.min(num_rows + k, num_cols)
    end

    select_diag = spzeros(T, sz_diag, length(x.children[1]))
    for i in 1:sz_diag
        select_diag[i, start_index] = 1
        start_index += num_rows + 1
    end

    child_obj = conic_form!(context, only(children(x)))
    obj = operate(add_operation, T, sign(x), select_diag, child_obj)
    return obj
end
