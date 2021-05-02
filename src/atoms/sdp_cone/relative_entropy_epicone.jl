#############################################################################
# relative_entropy_epicone.jl
# Returns τ constrained to
#   e' * X^{1/2} * logm(X^{1/2}*Y^{-1}*X^{1/2}) * X^{1/2} * e \preceq τ
#
# This function implements the semidefinite programming approximation given in
# the reference below.  Parameters m and k control the accuracy of this
# approximation: m is the number of quadrature nodes to use and k the number
# of square-roots to take. See reference for more details.
#
# All expressions and atoms are subtypes of AbstractExpr.
# Please read expressions.jl first.
#
#REFERENCE
#   Ported from CVXQUAD which is based on the paper: "Semidefinite
#   approximations of matrix logarithm" by Hamza Fawzi, James Saunderson and
#   Pablo A. Parrilo (arXiv:1705.00812)
#############################################################################

struct RelativeEntropyEpiCone <: AbstractExpr
    head::Symbol
    id_hash::UInt64
    children::Tuple{AbstractExpr, AbstractExpr}
    size::Tuple{Int, Int}
    m::Integer
    k::Integer
    e::AbstractMatrix

    function RelativeEntropyEpiCone(X::AbstractExpr, Y::AbstractExpr, m::Integer, k::Integer, e::AbstractArray)
        children = (X, Y)
        if size(X) != size(Y)
            throw(DimensionMismatch("X and Y must be the same size"))
        end
        n = size(X)[1]
        if size(X) != (n, n)
            throw(DimensionMismatch("X and Y must be square"))
        end
        if size(e) == (n,)
            e = reshape(e, (n, 1))
        end
        if ndims(e) != 2 || size(e)[1] != n
            throw(DimensionMismatch("e matrix must have n rows"))
        end
        return new(:relative_entropy_epicone, hash(children), children, (n, n), m, k, e)
    end
end

function sign(atom::RelativeEntropyEpiCone)
    X = atom.children[1]
    Y = atom.children[2]
    # Result is PSD, but it looks like Semidefinite is supposed to be NoSign.
    if sign(X) == ComplexSign() || sign(Y) == ComplexSign() || sign(Constant(atom.e)) == ComplexSign()
        return ComplexSign()
    else
        return NoSign()
    end
end

function monotonicity(atom::RelativeEntropyEpiCone)
    return (NoMonotonicity(), NoMonotonicity())
end

function curvature(atom::RelativeEntropyEpiCone)
    return ConvexVexity()
end

function evaluate(atom::RelativeEntropyEpiCone)
    X = evaluate(atom.children[1])
    Y = evaluate(atom.children[2])
    rX = sqrt(X)
    return rX * log(rX * Y^-1 * rX) * rX
end

relative_entropy_epicone(X::AbstractExpr,   Y::AbstractExpr,   m::Integer, k::Integer, e::AbstractArray = Matrix(1.0*I, size(X))) = RelativeEntropyEpiCone(X, Y, m, k, e)
relative_entropy_epicone(X::AbstractMatrix, Y::AbstractExpr,   m::Integer, k::Integer, e::AbstractArray = Matrix(1.0*I, size(X))) = RelativeEntropyEpiCone(Constant(X), Y, m, k, e)
relative_entropy_epicone(X::AbstractExpr,   Y::AbstractMatrix, m::Integer, k::Integer, e::AbstractArray = Matrix(1.0*I, size(X))) = RelativeEntropyEpiCone(X, Constant(Y), m, k, e)
relative_entropy_epicone(X::AbstractMatrix, Y::AbstractMatrix, m::Integer, k::Integer, e::AbstractArray = Matrix(1.0*I, size(X))) = RelativeEntropyEpiCone(Constant(X), Constant(Y), m, k, e)

function glquad(m)
    # Compute Gauss-Legendre quadrature nodes and weights on [0, 1]
    # Code below is from [Trefethen, "Is Gauss quadrature better than
    # Clenshaw-Curtis?", SIAM Review 2008] and computes the weights and
    # nodes on [-1, 1].
    beta = [ 0.5 ./ sqrt(1-(2*i)^-2) for i in 1:(m-1) ] # 3-term recurrence coeffs
    T = diagm(1 => beta, -1 => beta) # Jacobi matrix
    s, V = eigen(T)
    w = 2*V[1, :].^2 # weights
    # Translate and scale to [0, 1]
    s = (s .+ 1)/2
    w = w'/2
    return s, w
end

function conic_form!(atom::RelativeEntropyEpiCone, unique_conic_forms)
    if !has_conic_form(unique_conic_forms, atom)
        X = atom.children[1]
        Y = atom.children[2]
        m = atom.m
        k = atom.k
        e = atom.e
        n = size(X)[1]
        r = size(e)[2]

        s, w = glquad(m)

        if sign(atom) == ComplexSign()
            τ = ComplexVariable(r, r)
            Z = ComplexVariable(n, n)
            T = [ ComplexVariable(r, r) for i in 1:m ]
        else
            τ = Variable(r, r)
            Z = Variable(n, n)
            T = [ Variable(r, r) for i in 1:m ]
        end

        conic_form!(Z in GeomMeanHypoCone(X, Y, 1//(2^k), false), unique_conic_forms)

        for ii=1:m
            # Note that we are dividing by w here because it is easier
            # to do this than to do sum w_i T(:,...,:,ii) later (cf. line that
            # involves τ)

            conic_form!(
                [e'*X*e - s[ii]*T[ii]/w[ii]   e'*X;
                    X*e                       (1-s[ii])*X+s[ii]*Z] ⪰ 0,
                unique_conic_forms)

        end

        conic_form!((2^k)*sum(T) + τ ⪰ 0, unique_conic_forms)

        cache_conic_form!(unique_conic_forms, atom, conic_form!(τ, unique_conic_forms))
    end
    return get_conic_form(unique_conic_forms, atom)
end
