#############################################################################
# lieb_ando.jl
# Returns tr(K' * A^{1-t} * K * B^t) where A and B are positive semidefinite
# matrices and K is an arbitrary matrix (possibly rectangular).
#
# Disciplined convex programming information:
#    lieb_ando(A,B,K,t) is concave in (A,B) for t in [0,1], and convex
#    in (A,B) for t in [-1,0] or [1,2]. K is a fixed matrix.
#
# All expressions and atoms are subtypes of AbstractExpr.
# Please read expressions.jl first.
#
#REFERENCE
#   Ported from CVXQUAD which is based on the paper: "Lieb's concavity
#   theorem, matrix geometric means and semidefinite optimization" by Hamza
#   Fawzi and James Saunderson (arXiv:1512.03401)
#############################################################################

const MatrixOrConstant = Union{AbstractMatrix, Constant}

function lieb_ando(A::MatrixOrConstant, B::MatrixOrConstant, K::MatrixOrConstant, t::Rational)
    return real(tr(K' * A^(1-t) * K * B^t))
end

function lieb_ando(A::MatrixOrConstant, B::AbstractExpr, K::MatrixOrConstant, t::Rational)
    KAK = K' * A^(1-t) * K
    KAK = (KAK+KAK')/2
    return trace_mpower(B, t, KAK)
end

function lieb_ando(A::AbstractExpr, B::MatrixOrConstant, K::MatrixOrConstant, t::Rational)
    KBK = K * B^t * K'
    KBK = (KBK+KBK')/2
    return trace_mpower(A, 1-t, KBK)
end

function lieb_ando(A::AbstractExpr, B::AbstractExpr, K::MatrixOrConstant, t::Rational)
    n = size(A,1)
    m = size(B,1)
    Kvec = reshape(K',n*m,1)
    KvKv = Kvec * Kvec'
    KvKv = (KvKv+KvKv')/2
    Im = Matrix(1.0*I, m, m)
    In = Matrix(1.0*I, n, n)

    is_complex = sign(A) == ComplexSign() || sign(B) == ComplexSign() || sign(Constant(K)) == ComplexSign()
    if is_complex
        T = HermitianSemidefinite(n*m)
    else
        T = Semidefinite(n*m)
    end

    if t >= 0 && t <= 1
        # Concave function
        add_constraint!(T, T in GeomMeanHypoCone(kron(A,Im), kron(In,conj(B)), t, false))
        return real(tr(KvKv * T))
    elseif (t >= -1 && t <= 0) || (t >= 1 && t <= 2)
        # Convex function
        add_constraint!(T, T in GeomMeanEpiCone(kron(A,Im), kron(In,conj(B)), t, false))
        return real(tr(KvKv * T))
    else
        throw(DomainError(t, "t must be between -1 and 2"))
    end
end
