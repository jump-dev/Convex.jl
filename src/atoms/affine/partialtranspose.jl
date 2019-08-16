import Base.sign
export partialtranspose, PartialTransposeAtom
export sign, curvature, monotonicity, evaluate, conic_form!

struct PartialTransposeAtom <: AbstractExpr
    head::Symbol
    id_hash::UInt64
    children::Tuple{AbstractExpr}
    size::Tuple{Int, Int}
    sys::Int
    dims::Vector

    function PartialTransposeAtom(x::AbstractExpr, sys::Int, dims::Vector)
        if x.size[1] ≠ x.size[2]
            error("Only square matrices are supported")
        end
        if ! (1 ≤ sys ≤ length(dims))
            error("Invalid system, should between 1 and ", length(dims), "; got ", sys)
        end
        if x.size[1] ≠ prod(dims)
            error("Dimension of system doesn't correspond to dimension of subsystems")
        end
        children = (x, )
        return new(:partialtranspose, hash(children), children, x.size, sys, dims)    
    end
end

function sign(x::PartialTransposeAtom)
    return sign(x.children[1])
end

function curvature(x::PartialTransposeAtom)
    return ConstVexity()
end

function monotonicity(x::PartialTransposeAtom)
    return (Nondecreasing(),)
end

function evaluate(x::PartialTransposeAtom)
    return partialtranspose( evaluate( x.children[1] ))
end

"""
    partialtranspose(x, sys::Int, dims::Vector)

Returns the partial transpose of `x` over the `sys`th system, where `dims` is a vector of integers encoding the dimensions of each subsystem.
"""
function partialtranspose(x::AbstractMatrix, sys::Int, dims::Vector) 
    n = length(dims)
    d = prod(dims)
    s = n - sys + 1
    p = collect(1:2n)
    p[s] = n + s
    p[n + s] = s

    rdims = reverse(dims)
    r = reshape(x, [rdims;rdims]...)
    return reshape(permutedims(r,p),(d,d))
end


function permutedims_matrix(dims, p)
    d = prod(dims)
    n = length(dims)
    sparse(reshape(PermutedDimsArray(reshape(Matrix(I,d,d), (dims..., dims...)), (p..., (n + 1:2n)...)), (d, d)))
end
    
# borrowing this from transpose.jl: 
# Since everything is vectorized, we simply need to multiply x by a permutation
# matrix such that coeff * vectorized(x) - vectorized(x') = 0
function conic_form!(x::PartialTransposeAtom, unique_conic_forms::UniqueConicForms=UniqueConicForms())
    if !has_conic_form(unique_conic_forms, x)
        objective = conic_form!(x.children[1], unique_conic_forms)

        n = length(x.dims)
        d = prod(x.dims)
        s = n - x.sys + 1
        p = collect(1:2n)
        p[s] = n + s
        p[n + s] = s

        rdims = reverse(x.dims)

        partialtranspose_matrix = permutedims_matrix([rdims;rdims],p)

        objective = partialtranspose_matrix * objective
        cache_conic_form!(unique_conic_forms, x, objective)
    end
    return get_conic_form(unique_conic_forms, x)
end

partialtranspose(x::AbstractExpr,sys::Int, dim::Vector) = PartialTransposeAtom(x,sys,dim)


