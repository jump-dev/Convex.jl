#############################################################################
# trace_mpower(A, t, C) returns tr(C*A^t) where A and C are positive definite
# matrices and C is constant and t ∈ [-1, 2].
#
# When t ∈ [0,1], trace_mpower(A, t, C) is concave in A (for fixed
# positive semidefinite matrix C) and convex for t ∈ [-1,0] or [1,2].
#
# All expressions and atoms are subtypes of AbstractExpr.
# Please read expressions.jl first.
#
#REFERENCE
#   Ported from CVXQUAD which is based on the paper: "Lieb's concavity
#   theorem, matrix geometric means and semidefinite optimization" by Hamza
#   Fawzi and James Saunderson (arXiv:1512.03401)
#############################################################################

struct TraceMpower <: AbstractExpr
    head::Symbol
    id_hash::UInt64
    children::Tuple{AbstractExpr}
    size::Tuple{Int,Int}
    C::AbstractMatrix
    t::Rational

    function TraceMpower(A::AbstractExpr, t::Rational, C::AbstractMatrix)
        children = (A,)
        if size(A) != size(C)
            throw(DimensionMismatch("A and C must be the same size"))
        end
        n = size(A)[1]
        if size(A) != (n, n)
            throw(DimensionMismatch("A and C must be square"))
        end
        if norm(C - C') > 1e-6
            throw(DomainError(C, "C must be Hermitian"))
        end
        if any(eigvals(Hermitian(C)) .< -1e-6)
            throw(DomainError(C, "C must be positive semidefinite"))
        end
        if t < -1 || t > 2
            throw(DomainError(t, "t must be in the range [-1, 2]"))
        end
        return new(:trace_mpower, hash(children), children, (1, 1), C, t)
    end
end

function sign(atom::TraceMpower)
    return NoSign()
end

function monotonicity(atom::TraceMpower)
    return (NoMonotonicity(),)
end

function curvature(atom::TraceMpower)
    t = atom.t
    if t >= 0 && t <= 1
        return ConcaveVexity()
    else
        return ConvexVexity()
    end
end

function evaluate(atom::TraceMpower)
    A = evaluate(atom.children[1])
    return trace_mpower(A, atom.t, atom.C)
end

const MatrixOrConstant = Union{AbstractMatrix,Constant}

function trace_mpower(A::AbstractExpr, t::Rational, C::MatrixOrConstant)
    #println("trace_mpower general case")
    return TraceMpower(A, t, evaluate(C))
end

function trace_mpower(A::MatrixOrConstant, t::Rational, C::MatrixOrConstant)
    #println("trace_mpower constant A")
    return tr(C * A^t)
end

function trace_mpower(A::AbstractExprOrValue, t::Integer, C::MatrixOrConstant)
    return trace_mpower(A, t // 1, C)
end

function conic_form!(context::Context, atom::TraceMpower)
    A = atom.children[1]
    C = atom.C
    t = atom.t
    eye = Matrix(1.0 * I, size(A))

    is_complex = sign(A) == ComplexSign()
    if is_complex
        make_temporary = () -> HermitianSemidefinite(size(A)[1])
    else
        make_temporary = () -> Semidefinite(size(A)[1])
    end

    T = make_temporary()

    if t >= 0 && t <= 1
        add_constraints_to_context(
            T in GeomMeanHypoCone(eye, A, t, false),
            context,
        )
        # It's already a real mathematically, but need to make it a real type.
        u = real(tr(C * T))
        return conic_form!(context, maximize(u))
    else
        add_constraints_to_context(
            T in GeomMeanEpiCone(eye, A, t, false),
            context,
        )
        # It's already a real mathematically, but need to make it a real type.
        u = real(tr(C * T))
        return conic_form!(context, minimize(u))
    end
end
