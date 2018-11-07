#############################################################################
# lambdamin_max.jl
# Handles maximum and minimum eigenvalue of a symmetric positive definite matrix
# (and imposes the constraint that its argument be PSD)
# All expressions and atoms are subtypes of AbstractExpr.
# Please read expressions.jl first.
#############################################################################
export lambdamax, lambdamin

### Lambda max

struct LambdaMaxAtom <: AbstractExpr
    head::Symbol
    id_hash::UInt64
    children::Tuple{AbstractExpr}
    size::Tuple{Int, Int}

    function LambdaMaxAtom(x::AbstractExpr)
        children = (x,)
        m, n = size(x)
        if m == n
            return new(:lambdamax, hash(children), children, (1, 1))
        else
            error("lambdamax can only be applied to a square matrix.")
        end
    end
end

function sign(x::LambdaMaxAtom)
    return Positive()
end

function monotonicity(x::LambdaMaxAtom)
    return (Nondecreasing(),)
end

function curvature(x::LambdaMaxAtom)
    return ConvexVexity()
end

function evaluate(x::LambdaMaxAtom)
    eigvals(evaluate(x.children[1]))[end]
end

lambdamax(x::AbstractExpr) = LambdaMaxAtom(x)

# Create the equivalent conic problem:
#   minimize t
#   subject to
#            tI - A is positive semidefinite
#            A      is positive semidefinite
function conic_form!(x::LambdaMaxAtom, unique_conic_forms)
    if !has_conic_form(unique_conic_forms, x)
        A = x.children[1]
        m, n = size(A)
        t = Variable()
        p = minimize(t, t*Matrix(1.0I, n, n) - A ⪰ 0)
        cache_conic_form!(unique_conic_forms, x, p)
    end
    return get_conic_form(unique_conic_forms, x)
end

### Lambda min

struct LambdaMinAtom <: AbstractExpr
    head::Symbol
    id_hash::UInt64
    children::Tuple{AbstractExpr}
    size::Tuple{Int, Int}

    function LambdaMinAtom(x::AbstractExpr)
        children = (x,)
        m, n = size(x)
        if m == n
            return new(:lambdamin, hash(children), children, (1,1))
        else
            error("lambdamin can only be applied to a square matrix.")
        end
    end
end

function sign(x::LambdaMinAtom)
    return Positive()
end

function monotonicity(x::LambdaMinAtom)
    return (Nondecreasing(),)
end

function curvature(x::LambdaMinAtom)
    return ConcaveVexity()
end

function evaluate(x::LambdaMinAtom)
    eigvals(evaluate(x.children[1]))[1]
end

lambdamin(x::AbstractExpr) = LambdaMinAtom(x)

# Create the equivalent conic problem:
#   maximize t
#   subject to
#            A - tI is positive semidefinite
#            A      is positive semidefinite
function conic_form!(x::LambdaMinAtom, unique_conic_forms)
    if !has_conic_form(unique_conic_forms, x)
        A = x.children[1]
        m, n = size(A)
        t = Variable()
        p = maximize(t, A - t*Matrix(1.0I, n, n) ⪰ 0)
        cache_conic_form!(unique_conic_forms, x, p)
    end
    return get_conic_form(unique_conic_forms, x)
end
