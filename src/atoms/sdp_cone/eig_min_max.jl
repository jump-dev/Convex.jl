#############################################################################
# eigmin_max.jl
# Handles maximum and minimum eigenvalue of a symmetric positive definite matrix
# (and imposes the constraint that its argument be PSD)
# All expressions and atoms are subtypes of AbstractExpr.
# Please read expressions.jl first.
#############################################################################

### Eig max

mutable struct EigMaxAtom <: AbstractExpr
    children::Tuple{AbstractExpr}
    size::Tuple{Int,Int}

    function EigMaxAtom(x::AbstractExpr)
        children = (x,)
        m, n = size(x)
        if m == n
            return new(children, (1, 1))
        else
            error("eigmax can only be applied to a square matrix.")
        end
    end
end

head(io::IO, ::EigMaxAtom) = print(io, "eigmax")

function Base.sign(x::EigMaxAtom)
    return NoSign()
end

function monotonicity(x::EigMaxAtom)
    return (Nondecreasing(),)
end

function curvature(x::EigMaxAtom)
    return ConvexVexity()
end

function evaluate(x::EigMaxAtom)
    return LinearAlgebra.eigmax(evaluate(x.children[1]))
end

LinearAlgebra.eigmax(x::AbstractExpr) = EigMaxAtom(x)

# Create the equivalent conic problem:
#   minimize t
#   subject to
#            tI - A is positive semidefinite
#            A      is positive semidefinite
function new_conic_form!(context::Context{T}, x::EigMaxAtom) where {T}
    A = x.children[1]
    m, n = size(A)
    t = Variable()
    p = minimize(t, t * Matrix(one(T) * LinearAlgebra.I, n, n) - A ⪰ 0)
    return conic_form!(context, p)
end

### Eig min

mutable struct EigMinAtom <: AbstractExpr
    children::Tuple{AbstractExpr}
    size::Tuple{Int,Int}

    function EigMinAtom(x::AbstractExpr)
        children = (x,)
        m, n = size(x)
        if m == n
            return new(children, (1, 1))
        else
            error("eigmin can only be applied to a square matrix.")
        end
    end
end

head(io::IO, ::EigMinAtom) = print(io, "eigmin")

function Base.sign(x::EigMinAtom)
    return NoSign()
end

function monotonicity(x::EigMinAtom)
    return (Nondecreasing(),)
end

function curvature(x::EigMinAtom)
    return ConcaveVexity()
end

function evaluate(x::EigMinAtom)
    return LinearAlgebra.eigmin(evaluate(x.children[1]))
end

LinearAlgebra.eigmin(x::AbstractExpr) = EigMinAtom(x)

# Create the equivalent conic problem:
#   maximize t
#   subject to
#            A - tI is positive semidefinite
#            A      is positive semidefinite
function new_conic_form!(context::Context, x::EigMinAtom)
    A = x.children[1]
    m, n = size(A)
    t = Variable()
    p = maximize(t, A - t * Matrix(1.0 * LinearAlgebra.I, n, n) ⪰ 0)
    return conic_form!(context, p)
end
