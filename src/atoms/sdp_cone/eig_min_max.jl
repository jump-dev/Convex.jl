#############################################################################
# eigmin_max.jl
# Handles maximum and minimum eigenvalue of a symmetric positive definite matrix
# (and imposes the constraint that its argument be PSD)
# All expressions and atoms are subtypes of AbstractExpr.
# Please read expressions.jl first.
#############################################################################

import LinearAlgebra: eigmin, eigmax

### Eig max

struct EigMaxAtom <: AbstractExpr
    head::Symbol
    id_hash::UInt64
    children::Tuple{AbstractExpr}
    size::Tuple{Int,Int}

    function EigMaxAtom(x::AbstractExpr)
        children = (x,)
        m, n = size(x)
        if m == n
            return new(:eigmax, hash(children), children, (1, 1))
        else
            error("eigmax can only be applied to a square matrix.")
        end
    end
end

function sign(x::EigMaxAtom)
    return NoSign()
end

function monotonicity(x::EigMaxAtom)
    return (Nondecreasing(),)
end

function curvature(x::EigMaxAtom)
    return ConvexVexity()
end

function evaluate(x::EigMaxAtom)
    return eigmax(evaluate(x.children[1]))
end

eigmax(x::AbstractExpr) = EigMaxAtom(x)

# Create the equivalent conic problem:
#   minimize t
#   subject to
#            tI - A is positive semidefinite
#            A      is positive semidefinite
function conic_form!(x::EigMaxAtom, unique_conic_forms)
    if !has_conic_form(unique_conic_forms, x)
        A = x.children[1]
        m, n = size(A)
        t = Variable()
        p = minimize(t, t * Matrix(1.0I, n, n) - A ⪰ 0)
        cache_conic_form!(unique_conic_forms, x, p)
    end
    return get_conic_form(unique_conic_forms, x)
end

### Eig min

struct EigMinAtom <: AbstractExpr
    head::Symbol
    id_hash::UInt64
    children::Tuple{AbstractExpr}
    size::Tuple{Int,Int}

    function EigMinAtom(x::AbstractExpr)
        children = (x,)
        m, n = size(x)
        if m == n
            return new(:eigmin, hash(children), children, (1, 1))
        else
            error("eigmin can only be applied to a square matrix.")
        end
    end
end

function sign(x::EigMinAtom)
    return NoSign()
end

function monotonicity(x::EigMinAtom)
    return (Nondecreasing(),)
end

function curvature(x::EigMinAtom)
    return ConcaveVexity()
end

function evaluate(x::EigMinAtom)
    return eigmin(evaluate(x.children[1]))
end

eigmin(x::AbstractExpr) = EigMinAtom(x)

# Create the equivalent conic problem:
#   maximize t
#   subject to
#            A - tI is positive semidefinite
#            A      is positive semidefinite
function conic_form!(x::EigMinAtom, unique_conic_forms)
    if !has_conic_form(unique_conic_forms, x)
        A = x.children[1]
        m, n = size(A)
        t = Variable()
        p = maximize(t, A - t * Matrix(1.0I, n, n) ⪰ 0)
        cache_conic_form!(unique_conic_forms, x, p)
    end
    return get_conic_form(unique_conic_forms, x)
end
