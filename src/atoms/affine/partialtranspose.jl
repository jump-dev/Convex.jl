"""
    partialtranspose(x, sys::Int, dims::Vector)

Returns the partial transpose of `x` over the `sys`th system, where `dims` is a
vector of integers encoding the dimensions of each subsystem.
"""
function partialtranspose(x::AbstractMatrix, sys::Int, dims::Vector)
    if size(x, 1) != size(x, 2)
        throw(ArgumentError("Only square matrices are supported"))
    end
    if !(1 <= sys <= length(dims))
        msg = "Invalid system, should between 1 and $(length(dims)); got $sys"
        throw(ArgumentError(msg))
    end
    if size(x, 1) != prod(dims)
        msg = "Dimension of system doesn't correspond to dimension of subsystems"
        throw(ArgumentError(msg))
    end
    n = length(dims)
    d = prod(dims)
    s = n - sys + 1
    p = collect(1:2n)
    p[s] = n + s
    p[n+s] = s
    rdims = reverse(dims)
    r = reshape(x, (rdims..., rdims...))
    return reshape(permutedims(r, p), (d, d))
end

"""
    permutedims_matrix(dims, p)

Returns a matrix `M` so that for any vector `v` of length `prod(dims)`,
`M*v == vec(permutedims(reshape(v, dims), p))`.
"""
function permutedims_matrix(dims, p)
    d, n = prod(dims), length(dims)
    return SparseArrays.sparse(
        reshape(
            PermutedDimsArray(
                reshape(Matrix(LinearAlgebra.I, d, d), (dims..., dims...)),
                (p..., (n+1:2n)...),
            ),
            (d, d),
        ),
    )
end

# Vectorized implementation of the above
function partialtranspose(x::AbstractExpr, sys::Int, dims::Vector)
    if size(x, 1) != size(x, 2)
        throw(ArgumentError("Only square matrices are supported"))
    end
    if !(1 <= sys <= length(dims))
        msg = "Invalid system, should between 1 and $(length(dims)); got $sys"
        throw(ArgumentError(msg))
    end
    if size(x, 1) != prod(dims)
        msg = "Dimension of system doesn't correspond to dimension of subsystems"
        throw(ArgumentError(msg))
    end
    n = length(dims)
    d = prod(dims)
    s = n - sys + 1
    p = collect(1:2n)
    p[s] = n + s
    p[n+s] = s
    rdims = reverse(dims)
    partialtranspose_matrix = permutedims_matrix((rdims..., rdims...), p)
    return reshape(partialtranspose_matrix * vec(x), size(x)...)
end
