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

mutable struct Counts
    seen::Base.IdSet
    n_variables::Int
    n_scalar_variables::Int
    n_constraints::Int
    n_scalar_constraints::Int
    n_atoms::Int
    n_nonzeros::Int
    Counts() = new(Base.IdSet(), 0, 0, 0, 0, 0, 0)
end

function Counts(p::Problem)
    counts = Counts()
    _add_to_problem_count(counts, p)
    counts.n_atoms -= 1  # Don't count p as an atom
    return counts
end

function _add_to_problem_count(counts::Counts, node::AbstractExpr)
    if node in counts.seen
        return
    end
    for n in AbstractTrees.PreOrderDFS(node)
        if n !== node && !(n in counts.seen)
            _add_to_problem_count(counts, n)
        end
    end
    counts.n_atoms += 1
    push!(counts.seen, node)
    return
end

function _add_to_problem_count(counts::Counts, node::AbstractVariable)
    if node in counts.seen
        return
    end
    counts.n_variables += 1
    counts.n_scalar_variables += (iscomplex(node) ? 2 : 1) * length(node)
    push!(counts.seen, node)
    for c in get_constraints(node)
        _add_to_problem_count(counts, c)
    end
    return
end

_add_to_problem_count(::Counts, ::Vector{Constraint}) = nothing

_add_to_problem_count(::Counts, ::Nothing) = nothing

function _add_to_problem_count(counts::Counts, node::GenericConstraint)
    counts.n_constraints += 1
    counts.n_scalar_constraints +=
        (iscomplex(node) ? 2 : 1) * MOI.dimension(node.set)
    _add_to_problem_count(counts, node.child)
    return
end

function _add_to_problem_count(counts::Counts, node::Constraint)
    counts.n_constraints += 1
    # TODO(odow): we don't now the general case
    return
end

function _add_to_problem_count(counts::Counts, node::Constant)
    _add_to_problem_count(counts, node.value)
    push!(counts.seen, node)
    return
end

function _add_to_problem_count(counts::Counts, node::ComplexConstant)
    _add_to_problem_count(counts, node.real_constant)
    _add_to_problem_count(counts, node.imag_constant)
    push!(counts.seen, node)
    return
end

function _add_to_problem_count(counts::Counts, node::Number)
    counts.n_nonzeros += iscomplex(node) ? 2 : 1
    return
end

function _add_to_problem_count(counts::Counts, node::AbstractArray{<:Number})
    counts.n_nonzeros += length(node)
    return
end

function _add_to_problem_count(
    counts::Counts,
    node::SparseArrays.SparseMatrixCSC
)
    counts.n_nonzeros += SparseArrays.nnz(node)
    return
end

