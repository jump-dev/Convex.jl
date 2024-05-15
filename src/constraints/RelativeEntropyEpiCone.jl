# Copyright (c) 2014: Madeleine Udell and contributors
# Copyright (c) 2021: Hamza Fawzi
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

"""
    RelativeEntropyEpiConeSquare(
        side_dimension::Int,
        m::Integer = 3,
        k::Integer = 3,
        e::AbstractArray = Matrix(1.0 * LinearAlgebra.I(side_dimension)),
    )

Constrains `(τ, X, Y)` to:
```
τ ⪰ e' * X^{1/2} * logm(X^{1/2}*Y^{-1}*X^{1/2}) * X^{1/2} * e
```

This set implements the semidefinite programming approximation given in the
reference below.

Parameters `m` and `k` control the accuracy of this approximation: `m` is the
number of quadrature nodes to use and `k` the number of square-roots to take.
See reference for more details.

## Reference

Ported from CVXQUAD which is based on the paper: "Semidefinite approximations of
matrix logarithm" by Hamza Fawzi, James Saunderson and Pablo A. Parrilo
(arXiv:1705.00812)
"""
struct RelativeEntropyEpiConeSquare <: MOI.AbstractVectorSet
    side_dimension::Int
    m::Integer
    k::Integer
    e::AbstractMatrix

    function RelativeEntropyEpiConeSquare(
        side_dimension::Int,
        m::Integer = 3,
        k::Integer = 3,
        e::AbstractArray = Matrix(1.0 * LinearAlgebra.I(side_dimension)),
    )
        if length(size(e)) == 1
            e = reshape(e, (size(e, 1), 1))
        end
        if ndims(e) != 2 || size(e, 1) != n
            throw(DimensionMismatch("e matrix must have n rows"))
        end
        return new(side_dimension, m, k, e)
    end
end

function MOI.dimension(set::RelativeEntropyEpiConeSquare)
    return size(set.e, 2)^2 + 2 * set.side_dimension^2
end

function head(io::IO, ::RelativeEntropyEpiConeSquare)
    return print(io, "RelativeEntropyEpiCone")
end

function _get_matrices(c::GenericConstraint{RelativeEntropyEpiConeSquare})
    n_τ, n_x = size(c.set.e, 2), c.set.side_dimension
    d_τ, d_x = n_τ^2, n_x^2
    τ = reshape(c.child[1:d_τ], n_τ, n_τ)
    X = reshape(c.child[d_τ.+(1:d_x)], n_x, n_x)
    Y = reshape(c.child[d_τ.+d_x.+(1:d_x)], n_x, n_x)
    return τ, X, Y
end

# This negative relative entropy function is matrix convex (arxiv:1705.00812).
# So if X and Y are convex sets, then τ ⪰ -D_op(X || Y) will be a convex set.
function vexity(constraint::GenericConstraint{RelativeEntropyEpiConeSquare})
    τ, X, Y = _get_matrices(constraint)
    if vexity(X) in (ConcaveVexity(), NotDcp()) ||
       vexity(Y) in (ConcaveVexity(), NotDcp())
        return NotDcp()
    end
    return ConvexVexity() + -vexity(τ)
end

"""
Compute Gauss-Legendre quadrature nodes and weights on [0, 1].

Code below is from Trefethen (2008), "Is Gauss quadrature better than
Clenshaw-Curtis?", SIAM Review, and computes the weights and nodes on [-1, 1].
"""
function _gauss_legendre_quadrature(m)
    # 3-term recurrence coeffs
    beta = [0.5 ./ sqrt(1 - (2 * i)^-2) for i in 1:(m-1)]
    # Jacobi matrix
    T = LinearAlgebra.diagm(1 => beta, -1 => beta)
    s, V = LinearAlgebra.eigen(T)
    w = 2 * V[1, :] .^ 2
    # Translate and scale to [0, 1]
    s = (s .+ 1) / 2
    w = w' / 2
    return s, w
end

function _add_constraint!(
    context::Context,
    constraint::GenericConstraint{RelativeEntropyEpiConeSquare},
)
    τ, X, Y = _get_matrices(constraint)
    m, k, e = constraint.set.m, constraint.set.k, constraint.set.e
    n = size(X, 1)
    r = size(e, 2)
    s, w = _gauss_legendre_quadrature(m)
    is_complex =
        sign(X) == ComplexSign() ||
        sign(Y) == ComplexSign() ||
        sign(constant(e)) == ComplexSign()
    if is_complex
        Z = ComplexVariable(n, n)
        T = [ComplexVariable(r, r) for i in 1:m]
    else
        Z = Variable(n, n)
        T = [Variable(r, r) for i in 1:m]
    end
    add_constraint!(
        context,
        GenericConstraint(
            (Z, X, Y),
            GeometricMeanHypoConeSquare(1 // (2^k), n, false),
        ),
    )
    for ii in 1:m
        # Note that we are dividing by w here because it is easier to do this
        # than to do sum w_i T(:,...,:,ii) later (cf. line that involves τ)
        A_ii = [
            (e' * X * e - s[ii] * T[ii] / w[ii]) (e' * X)
            (X * e) ((1 - s[ii]) * X + s[ii] * Z)
        ]
        add_constraint!(context, A_ii ⪰ 0)
    end
    add_constraint!(context, 2^k * sum(T) + τ ⪰ 0)
    return
end
