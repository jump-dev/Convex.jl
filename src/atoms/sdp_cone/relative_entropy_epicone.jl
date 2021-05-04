#############################################################################
# relative_entropy_epicone.jl
# Constrains τ to
#   τ ⪰ e' * X^{1/2} * logm(X^{1/2}*Y^{-1}*X^{1/2}) * X^{1/2} * e
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

struct RelativeEntropyEpiCone
    X::AbstractExpr
    Y::AbstractExpr
    m::Integer
    k::Integer
    e::AbstractMatrix
    size::Tuple{Int, Int}

    function RelativeEntropyEpiCone(X::AbstractExpr, Y::AbstractExpr, m::Integer=3, k::Integer=3, e::AbstractArray = Matrix(1.0*I, size(X)))
        if size(X) != size(Y)
            throw(DimensionMismatch("X and Y must be the same size"))
        end
        n = size(X)[1]
        if size(X) != (n, n)
            throw(DimensionMismatch("X and Y must be square"))
        end
        if length(size(e)) == 1
            e = reshape(e, (size(e)[1], 1))
        end
        erows, ecols = size(e)
        if ndims(e) != 2 || erows != n
            throw(DimensionMismatch("e matrix must have n rows"))
        end
        return new(X, Y, m, k, e, (ecols, ecols))
    end

    RelativeEntropyEpiCone(X::Value,        Y::AbstractExpr, m::Integer=3, k::Integer=3, e::AbstractArray = Matrix(1.0*I, size(X))) = RelativeEntropyEpiCone(Constant(X), Y, m, k, e)
    RelativeEntropyEpiCone(X::AbstractExpr, Y::Value,        m::Integer=3, k::Integer=3, e::AbstractArray = Matrix(1.0*I, size(X))) = RelativeEntropyEpiCone(X, Constant(Y), m, k, e)
    RelativeEntropyEpiCone(X::Value,        Y::Value,        m::Integer=3, k::Integer=3, e::AbstractArray = Matrix(1.0*I, size(X))) = RelativeEntropyEpiCone(Constant(X), Constant(Y), m, k, e)
end

struct RelativeEntropyEpiConeConstraint <: Constraint
    head::Symbol
    id_hash::UInt64
    τ::AbstractExpr
    cone::RelativeEntropyEpiCone

    function RelativeEntropyEpiConeConstraint(τ::AbstractExpr, cone::RelativeEntropyEpiCone)
        if size(τ) != cone.size
            throw(DimensionMismatch("τ must be size $(cone.size)"))
        end
        id_hash = hash((cone.X, cone.Y, cone.m, cone.k, cone.e, :RelativeEntropyEpiCone))
        return new(:RelativeEntropyEpiCone, id_hash, τ, cone)
    end

    RelativeEntropyEpiConeConstraint(τ::Value, cone::RelativeEntropyEpiCone) = RelativeEntropyEpiConeConstraint(Constant(τ), cone)
end

in(τ, cone::RelativeEntropyEpiCone) = RelativeEntropyEpiConeConstraint(τ, cone)

function AbstractTrees.children(constraint::RelativeEntropyEpiConeConstraint)
    return (constraint.τ, constraint.cone.X, constraint.cone.Y)
end

# This negative relative entropy function is matrix convex (arxiv:1705.00812).
# So if X and Y are convex sets, then τ ⪰ -D_op(X || Y) will be a convex set.
function vexity(constraint::RelativeEntropyEpiConeConstraint)
    X = vexity(constraint.cone.X)
    Y = vexity(constraint.cone.Y)
    τ = vexity(constraint.τ)

    # NOTE: can't say X == NotDcp() because the NotDcp constructor prints a warning message.
    if typeof(X) == ConcaveVexity || typeof(X) == NotDcp
        return NotDcp()
    end
    if typeof(Y) == ConcaveVexity || typeof(Y) == NotDcp
        return NotDcp()
    end
    # Copied from vexity(c::GtConstraint)
    vex = ConvexVexity() + (-τ)
    if vex == ConcaveVexity()
        vex = NotDcp()
    end
    return vex
end

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

function conic_form!(constraint::RelativeEntropyEpiConeConstraint, unique_conic_forms::UniqueConicForms)
    if !has_conic_form(unique_conic_forms, constraint)
        X = constraint.cone.X
        Y = constraint.cone.Y
        m = constraint.cone.m
        k = constraint.cone.k
        e = constraint.cone.e
        τ = constraint.τ
        n = size(X)[1]
        r = size(e)[2]

        s, w = glquad(m)

        is_complex = sign(X) == ComplexSign() || sign(Y) == ComplexSign() || sign(Constant(e)) == ComplexSign()
        if is_complex
            Z = ComplexVariable(n, n)
            T = [ ComplexVariable(r, r) for i in 1:m ]
        else
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

        cache_conic_form!(unique_conic_forms, constraint, Array{Convex.ConicConstr,1}())
    end
    return get_conic_form(unique_conic_forms, constraint)
end
