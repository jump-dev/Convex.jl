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
