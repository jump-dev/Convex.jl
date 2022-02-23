#############################################################################
# transpose.jl
# Returns the transpose of a matrix
# All expressions and atoms are subtpyes of AbstractExpr.
# Please read expressions.jl first.
#############################################################################

import Base.transpose, Base.adjoint

struct TransposeAtom <: AbstractExpr
    head::Symbol
    id_hash::UInt64
    children::Tuple{AbstractExpr}
    size::Tuple{Int,Int}

    function TransposeAtom(x::AbstractExpr)
        children = (x,)
        return new(:transpose, hash(children), children, (x.size[2], x.size[1]))
    end
end

function sign(x::TransposeAtom)
    return sign(x.children[1])
end

function monotonicity(x::TransposeAtom)
    return (Nondecreasing(),)
end

function curvature(x::TransposeAtom)
    return ConstVexity()
end

function evaluate(x::TransposeAtom)
    return transpose(evaluate(x.children[1]))
end

# Since everything is vectorized, we simply need to multiply x by a permutation
# matrix such that coeff * vectorized(x) - vectorized(x') = 0
function conic_form!(x::TransposeAtom, unique_conic_forms::UniqueConicForms)
    if !has_conic_form(unique_conic_forms, x)
        objective = conic_form!(x.children[1], unique_conic_forms)

        sz = length(x)

        num_rows = x.size[1]
        num_cols = x.size[2]

        I = Array{Int}(undef, sz)
        J = Array{Int}(undef, sz)

        k = 1
        for r in 1:num_rows
            for c in 1:num_cols
                I[k] = (c - 1) * num_rows + r
                J[k] = (r - 1) * num_cols + c
                k += 1
            end
        end

        transpose_matrix = sparse(I, J, 1.0)

        objective = transpose_matrix * objective
        cache_conic_form!(unique_conic_forms, x, objective)
    end
    return get_conic_form(unique_conic_forms, x)
end

transpose(x::AbstractExpr) = TransposeAtom(x)

struct AdjointAtom <: AbstractExpr
    head::Symbol
    id_hash::UInt64
    children::Tuple{AbstractExpr}
    size::Tuple{Int,Int}

    function AdjointAtom(x::AbstractExpr)
        children = (x,)
        return new(:adjoint, hash(children), children, (x.size[2], x.size[1]))
    end
end

function sign(x::AdjointAtom)
    return sign(x.children[1])
end

function monotonicity(x::AdjointAtom)
    return (Nondecreasing(),)
end

function curvature(x::AdjointAtom)
    return ConstVexity()
end

function evaluate(x::AdjointAtom)
    return evaluate(x.children[1])'
end

# Since everything is vectorized, we simply need to multiply x by a permutation
# matrix such that coeff * vectorized(x) - vectorized(x') = 0
function conic_form!(x::AdjointAtom, unique_conic_forms::UniqueConicForms)
    if !has_conic_form(unique_conic_forms, x)
        objective = conic_form!(x.children[1], unique_conic_forms)

        sz = length(x)

        num_rows = x.size[1]
        num_cols = x.size[2]

        I = Array{Int}(undef, sz)
        J = Array{Int}(undef, sz)

        k = 1
        for r in 1:num_rows
            for c in 1:num_cols
                I[k] = (c - 1) * num_rows + r
                J[k] = (r - 1) * num_cols + c
                k += 1
            end
        end

        transpose_matrix = sparse(I, J, 1.0)
        objective = transpose_matrix * objective

        for var in keys(objective)
            x1 = conj(objective[var][1])
            x2 = conj(objective[var][2])
            objective[var] = (x1, x2)
        end
        cache_conic_form!(unique_conic_forms, x, objective)
    end
    return get_conic_form(unique_conic_forms, x)
end

adjoint(x::AbstractExpr) = AdjointAtom(x)
adjoint(x::Constant) = Constant(x.value')
