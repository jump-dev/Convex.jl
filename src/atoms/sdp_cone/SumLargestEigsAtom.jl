mutable struct SumLargestEigsAtom <: AbstractExpr
    children::Tuple{AbstractExpr,AbstractExpr}
    size::Tuple{Int,Int}

    function SumLargestEigsAtom(x::AbstractExpr, k::AbstractExpr)
        children = (x, k)
        m, n = size(x)
        if m != n
            error("sumlargesteigs can only be applied to a square matrix.")
        end
        return new(children, (1, 1))
    end
end

head(io::IO, ::SumLargestEigsAtom) = print(io, "sumlargesteigs")

Base.sign(::SumLargestEigsAtom) = NoSign()

monotonicity(::SumLargestEigsAtom) = (Nondecreasing(), NoMonotonicity())

curvature(::SumLargestEigsAtom) = ConvexVexity()

function evaluate(x::SumLargestEigsAtom)
    return LinearAlgebra.eigvals(evaluate(x.children[1]))[end-x.children[2]:end]
end

function sumlargesteigs(x::AbstractExpr, k::Int)
    return k == 0 ? Constant(0) : SumLargestEigsAtom(x, Constant(k))
end

# Create the equivalent conic problem:
#   minimize sk + Tr(Z)
#   subject to
#            Z ⪰ 0
#            A ⪰ 0
#            Z + sI ⪰ A
# See Ben-Tal and Nemirovski, "Lectures on Modern Convex Optimization"
# Example 18.c
function new_conic_form!(context::Context{T}, x::SumLargestEigsAtom) where {T}
    X, k = x.children
    m, n = size(X)
    Z = if iscomplex(sign(X))
        ComplexVariable(n, n)
    else
        Variable(n, n)
    end
    s = Variable()
    # Note: we know the trace is real, since Z is PSD, but we need to tell
    # Convex.jl that.
    p = minimize(
        s * k + real(LinearAlgebra.tr(Z)),
        [Z - X + s * Matrix(1.0 * LinearAlgebra.I, n, n) ⪰ 0, Z ⪰ 0, X == X'],
    )
    return conic_form!(context, p)
end
