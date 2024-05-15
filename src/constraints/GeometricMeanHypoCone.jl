# Copyright (c) 2014: Madeleine Udell and contributors
# Copyright (c) 2021: Hamza Fawzi
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

"""
Constrains T to
  A #_t B ⪰ T
where:
  * A #_t B is the t-weighted geometric mean of A and B:
    A^{1/2} (A^{-1/2} B A^{-1/2})^t A^{1/2}
    Parameter t should be in [0,1].
  * Constraints A ⪰ 0, B ⪰ 0 are added.

Note on parameter fullhyp:
  In many applications one doesn't need the full hypograph
    hyp_t = {(A,B,T) : A #_t B ⪰ T}
  but rather it is enough to work with a convex set C_t that satisfies
    (A,B,A #_t B) \\in C_t
    (A,B,T) \\in C_t  =>  A #_t B ⪰ T
  In this case one should set fullhyp = false. The SDP description will be
  (slightly) smaller. (By default fullhyp is set to true).

All expressions and atoms are subtypes of AbstractExpr.
Please read expressions.jl first.

REFERENCE
  Ported from CVXQUAD which is based on the paper: "Lieb's concavity
  theorem, matrix geometric means and semidefinite optimization" by Hamza
  Fawzi and James Saunderson (arXiv:1512.03401)
"""
struct GeometricMeanHypoConeSquare <: MOI.AbstractVectorSet
    t::Rational
    side_dimension::Int
    fullhyp::Bool

    function GeometricMeanHypoConeSquare(
        t::Rational,
        side_dimension::Int,
        fullhyp::Bool = true,
    )
        if t < 0 || t > 1
            throw(DomainError(t, "t must be in the range [0, 1]"))
        end
        return new(t, side_dimension, fullhyp)
    end
end

MOI.dimension(set::GeometricMeanHypoConeSquare) = 3 * set.side_dimension^2

function head(io::IO, ::GeometricMeanHypoConeSquare)
    return print(io, "GeometricMeanHypoConeSquare")
end

function Constraint(func::Tuple, set::GeometricMeanHypoConeSquare)
    for f in func
        n = LinearAlgebra.checksquare(f)
        if n != set.side_dimension
            throw(
                DimensionMismatch(
                    "Matrix of side dimension `$n` does not match set of side dimension `$(set.side_dimension)`",
                ),
            )
        end
    end
    return Constraint(vcat(vec.(func)...), set)
end

function _get_matrices(c::Constraint{GeometricMeanHypoConeSquare})
    n = c.set.side_dimension
    d = n^2
    T = reshape(c.child[1:d], n, n)
    A = reshape(c.child[d.+(1:d)], n, n)
    B = reshape(c.child[2d.+(1:d)], n, n)
    return T, A, B
end

# For t ∈ [0,1] the t-weighted matrix geometric mean is matrix concave (arxiv:1512.03401).
# So if A and B are convex sets, then T ⪯ A #_t B will be a convex set.
function vexity(constraint::Constraint{GeometricMeanHypoConeSquare})
    T, A, B = _get_matrices(constraint)
    if vexity(A) in (ConvexVexity(), NotDcp()) ||
       vexity(B) in (ConvexVexity(), NotDcp())
        return NotDcp()
    end
    return -ConcaveVexity() + vexity(T)
end

function _add_constraint!(
    context::Context,
    constraint::Constraint{GeometricMeanHypoConeSquare},
)
    T, A, B = _get_matrices(constraint)
    t = constraint.set.t
    fullhyp = constraint.set.fullhyp

    is_complex =
        sign(A) == ComplexSign() ||
        sign(B) == ComplexSign() ||
        sign(T) == ComplexSign()
    n = size(A, 1)
    if is_complex
        make_temporary = () -> HermitianSemidefinite(n)
    else
        make_temporary = () -> Semidefinite(n)
    end

    if fullhyp && t != 0 && t != 1
        W = make_temporary()
        add_constraint!(
            context,
            Constraint(
                (W, A, B),
                GeometricMeanHypoConeSquare(t, n, false),
            ),
        )
        add_constraint!(context, W ⪰ T)
    else
        p = t.num
        q = t.den

        if t == 0
            #println("geom_mean_hypocone p=$p q=$q t==0")
            add_constraint!(context, A ⪰ 0)
            add_constraint!(context, B ⪰ 0)
            add_constraint!(context, A ⪰ T)
        elseif t == 1
            #println("geom_mean_hypocone p=$p q=$q t==1")
            add_constraint!(context, A ⪰ 0)
            add_constraint!(context, B ⪰ 0)
            add_constraint!(context, B ⪰ T)
        elseif t == 1 // 2
            #println("geom_mean_hypocone p=$p q=$q t==1/2")
            add_constraint!(context, [A T; T' B] ⪰ 0)
        elseif ispow2(q)
            #println("geom_mean_hypocone p=$p q=$q ispow2(q)")
            Z = make_temporary()
            if t < 1 / 2
                add_constraint!(
                    context,
                    Constraint(
                        (Z, A, B),
                        GeometricMeanHypoConeSquare(2 * t, n, false),
                    ),
                )
                add_constraint!(context, [A T; T' Z] ⪰ 0)
            else
                add_constraint!(
                    context,
                    Constraint(
                        (Z, A, B),
                        GeometricMeanHypoConeSquare(2 * t - 1, n, false),
                    ),
                )
                add_constraint!(context, [B T; T' Z] ⪰ 0)
            end
        elseif ispow2(p) && t > 1 // 2
            #println("geom_mean_hypocone p=$p q=$q ispow2(p) && t>1/2")
            Z = make_temporary()
            add_constraint!(
                context,
                Constraint(
                    (Z, A, T),
                    GeometricMeanHypoConeSquare((2 * p - q) // p, n, false),
                ),
            )
            add_constraint!(context, [Z T; T B] ⪰ 0)
        elseif t < 1 / 2
            #println("geom_mean_hypocone p=$p q=$q t<1/2")
            X = make_temporary()
            # Decompose t = (p/2^l) * (2^l/q) where l=floor(log2(q))
            l = floor(Int, log2(q))
            add_constraint!(
                context,
                Constraint(
                    (X, A, B),
                    GeometricMeanHypoConeSquare(p // (2^l), n, false),
                ),
            )
            add_constraint!(
                context,
                Constraint(
                    (T, A, X),
                    GeometricMeanHypoConeSquare((2^l) // q, n, false),
                ),
            )
        else
            #println("geom_mean_hypocone p=$p q=$q else")
            add_constraint!(
                context,
                Constraint(
                    (T, B, A),
                    GeometricMeanHypoConeSquare(1 - t, n, false),
                ),
            )
        end
    end
    return
end
