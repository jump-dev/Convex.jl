export partialtrace

# We compute the partial trace of x by summing over
# (I ⊗ <j| ⊗ I) x (I ⊗ |j> ⊗ I) for all j's
# in the system we want to trace out.
# This function returns the jth term in the sum, namely
# (I ⊗ <j| ⊗ I) x (I ⊗ |j> ⊗ I).
function _term(x, j::Int, sys, dims)
    a = sparse(1.0I, 1, 1)
    b = sparse(1.0I, 1, 1)
    for (i_sys, dim) in enumerate(dims)
        if i_sys == sys
            # create a vector that is only 1 at its jth component
            v = spzeros(dim, 1);
            v[j] = 1;
            a = kron(a, v')
            b = kron(b, v)
        else
            a = kron(a, sparse(1.0I, dim, dim))
            b = kron(b, sparse(1.0I, dim, dim))
        end
    end
    return a * x * b
end

"""
    partialtrace(x, sys::Int, dims::Vector)

Returns the partial trace of `x` over the `sys`th system, where `dims` is a vector of integers encoding the dimensions of each subsystem.
"""
function partialtrace(x, sys::Int, dims::Vector)
    if size(x, 1) ≠ size(x, 2)
        throw(ArgumentError("only square matrices are supported"))
    end
    if ! (1 ≤ sys ≤ length(dims))
        throw(ArgumentError("invalid system index, should between 1 and $(length(dims)), got $sys"))
    end
    if size(x, 1) ≠ prod(dims)
        throw(ArgumentError("dimension of system doesn't correspond to dimension of subsystems"))
    end

    return sum(j -> _term(x, j, sys, dims), 1:dims[sys])
end