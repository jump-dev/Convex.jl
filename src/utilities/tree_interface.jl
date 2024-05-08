# Copyright (c) 2014: Madeleine Udell and contributors
#
# Use of this source code is governed by a BSD-style license that can be found
# in the LICENSE file or at https://opensource.org/license/bsd-2-clause

AbstractTrees.children(p::Problem) = (p.objective, p.constraints)

AbstractTrees.children(e::AbstractExpr) = e.children

AbstractTrees.children(v::AbstractVariable) = ()

AbstractTrees.children(c::Constant) = ()

AbstractTrees.children(C::Constraint) = (C.lhs, C.rhs)

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
                ")",
            )
        end
        print(io, ", ", underscorise(c.n_atoms), " atom", maybe_s(c.n_atoms))
    end
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

function counts_recursive(node)
    seen = Base.IdSet()
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
    return Counts(;
        n_constraints = 1,
        total_len_constraints = MOI.dimension(node.set),
    )
end

function counts(node::Union{RelativeEntropyEpiCone,GeometricMeanHypoCone})
    return Counts(; n_constraints = 1, total_len_constraints = node.size[1])
end

function counts(node::Union{Constant,ComplexConstant})
    f = iscomplex(node) ? 2 : 1
    return Counts(;
        n_constants = 1,
        total_len_dense_constants = SparseArrays.issparse(node.value) ? 0 :
                                    f * length(node.value),
        total_nnz_sparse_constants = SparseArrays.issparse(node.value) ?
                                     f * SparseArrays.nnz(node.value) : 0,
    )
end

# This node has no counts itself (we still count the children!)
counts(::Vector{Constraint}) = zero(Counts)
