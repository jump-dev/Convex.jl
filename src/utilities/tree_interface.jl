# Copyright (c) 2014: Madeleine Udell and contributors
#
# Use of this source code is governed by a BSD-style license that can be found
# in the LICENSE file or at https://opensource.org/license/bsd-2-clause

AbstractTrees.children(p::Problem) = (p.objective, p.constraints)

AbstractTrees.children(e::AbstractExpr) = e.children

AbstractTrees.children(v::AbstractVariable) = ()

AbstractTrees.children(c::Constant) = ()

AbstractTrees.children(C::Constraint) = (C.lhs, C.rhs)
AbstractTrees.children(C::RelativeEntropyEpiConeConstraint) = (C.τ, C.cone)

AbstractTrees.printnode(io::IO, node::AbstractExpr) = summary(io, node)

AbstractTrees.printnode(io::IO, node::Constraint) = summary(io, node)

function AbstractTrees.printnode(io::IO, node::Vector{<:Constraint})
    if length(node) == 0
        print(io, "no constraints")
    else
        print(io, "constraints")
    end
    return
end

AbstractTrees.printnode(io::IO, node::Problem) = print(io, node.head)

function AbstractTrees.printnode(io::IO, node::Constant)
    if length(node.value) <= 3
        show(IOContext(io, :compact => true), node.value)
    else
        summary(io, node.value)
    end
    return
end

AbstractTrees.printnode(io::IO, node::AbstractVariable) = summary(io, node)

Base.@kwdef struct Counts
    n_atoms::Int = 0
    n_variables::Int = 0
    n_constants::Int = 0
    n_constraints::Int = 0
    total_len_variables::Int = 0
    total_len_constraints::Int = 0
    total_len_dense_constants::Int = 0
    total_nnz_sparse_constants::Int = 0
end

# Borrowed from Flux
function underscorise(n::Integer)
    return join(
        reverse(join.(reverse.(Iterators.partition(digits(n), 3)))),
        '_',
    )
end

function Base.show(io::IO, c::Counts)
    maybe_s(n) = n == 1 ? "" : "s"
    print(io, c.n_variables, " variable", maybe_s(c.n_variables))
    if c.n_variables > 0
        print(
            io,
            " (",
            underscorise(c.total_len_variables),
            " element",
            maybe_s(c.total_len_variables),
            ")",
        )
    end

    print(
        io,
        ", ",
        underscorise(c.n_constraints),
        " constraint",
        maybe_s(c.n_constraints),
    )
    if c.n_constraints > 0
        print(
            io,
            " (",
            underscorise(c.total_len_constraints),
            " element",
            maybe_s(c.total_len_constraints),
            ")",
        )
    end
    print(
        io,
        ", ",
        underscorise(c.n_constants),
        " constant",
        maybe_s(c.n_constants),
    )
    if c.n_constants > 0
        print(io, " (")
        if c.total_len_dense_constants > 0
            print(
                io,
                underscorise(c.total_len_dense_constants),
                " dense element",
                maybe_s(c.total_len_dense_constants),
            )
        end
        if c.total_len_dense_constants > 0 && c.total_nnz_sparse_constants > 0
            print(io, ", ")
        end
        if c.total_nnz_sparse_constants > 0
            print(
                io,
                underscorise(c.total_nnz_sparse_constants),
                " sparse nonzero",
                maybe_s(c.total_nnz_sparse_constants),
            )
        end
        print(io, ")")
    end
    print(io, ", and ", underscorise(c.n_atoms), " atom", maybe_s(c.n_atoms))
    return nothing
end

function Base.:(+)(c1::Counts, c2::Counts)
    return Counts(
        (getfield(c1, i) + getfield(c2, i) for i in 1:fieldcount(Counts))...,
    )
end

Base.zero(::Type{Counts}) = Counts()

# Don't double-count
function single_count!(seen::Base.IdSet, node)
    if node in seen
        return zero(Counts)
    else
        push!(seen, node)
        return counts(node)
    end
end

# This is tricky: variables can carry their own constraints which are accessed via `get_constraints`.
# These constraints are *not* their children; often the constraint will involve
# the variable itself, so there would be a cycle (which would break the tree-based printing).
# However, we *do* want to count these constraints. Therefore, we will
# specially recurse on them here, whenever we encounter an `AbstractVariable`.
function single_count!(seen::Base.IdSet, node::AbstractVariable)
    if node in seen
        return zero(Counts)
    else
        push!(seen, node)
        return counts(node) + counts_recursive!(seen, get_constraints(node))
    end
end

function counts_recursive(node)
    return counts_recursive!(Base.IdSet(), node)
end

function counts_recursive!(seen, node)
    return sum(
        c -> single_count!(seen, c),
        AbstractTrees.PreOrderDFS(node);
        init = zero(Counts),
    )
end

# Here we define `counts` for each object.

function counts(::AbstractExpr)
    return Counts(; n_atoms = 1)
end

function counts(node::AbstractVariable)
    f = iscomplex(node) ? 2 : 1
    return Counts(; n_variables = 1, total_len_variables = f * length(node))
end

function counts(node::GenericConstraint)
    f = iscomplex(node) ? 2 : 1
    return Counts(;
        n_constraints = 1,
        total_len_constraints = f * MOI.dimension(node.set),
    )
end

function counts(node::GeometricMeanHypoCone)
    return Counts()
end

function counts(node::RelativeEntropyEpiCone)
    return Counts()
end

function counts(node::RelativeEntropyEpiConeConstraint)
    f = iscomplex(node) ? 2 : 1
    return Counts(; n_constraints = 1, total_len_constraints = f * node.cone.size[1])
end


function counts(node::Convex.GeometricMeanHypoConeConstraint)
    f = iscomplex(node) ? 2 : 1
    return Counts(; n_constraints = 1, total_len_constraints = f * node.cone.size[1])
end


# e.g. `satisfy` objectives
counts(node::Nothing) = Counts()

function counts(node::Constant)
    return Counts(;
        n_constants = 1,
        total_len_dense_constants = SparseArrays.issparse(node.value) ? 0 :
                                    length(node.value),
        total_nnz_sparse_constants = SparseArrays.issparse(node.value) ?
                                     SparseArrays.nnz(node.value) : 0,
    )
end

function counts(node::ComplexConstant)
    re = counts(node.real_constant)
    im = counts(node.imag_constant)
    return Counts(;
        n_constants = 1,
        total_len_dense_constants = re.total_len_dense_constants +
                                    im.total_len_dense_constants,
        total_nnz_sparse_constants = re.total_nnz_sparse_constants +
                                     im.total_nnz_sparse_constants,
    )
end

# This node has no counts itself (we still count the children!)
counts(::Vector{Constraint}) = zero(Counts)
