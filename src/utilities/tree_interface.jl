import AbstractTrees

function AbstractTrees.children(p::Problem)
    return (p.objective, p.constraints)
end

function AbstractTrees.children(e::AbstractExpr)
    return e.children
end

AbstractTrees.children(v::Variable) = ()
AbstractTrees.children(c::Constant) = ()

AbstractTrees.children(C::Constraint) = (C.lhs, C.rhs)

AbstractTrees.printnode(io::IO, node::AbstractExpr) = print(io,node.head)
AbstractTrees.printnode(io::IO, node::Constraint) = print(io,node.head)

AbstractTrees.printnode(io::IO, node::Vector{Constraint}) = print(io, "constraints")


AbstractTrees.printnode(io::IO, node::Problem) = print(io,node.head)

AbstractTrees.printnode(io::IO, node::Constant) = show(IOContext(io, :compact => true), node.value)

AbstractTrees.printnode(io::IO, node::Variable) = summary(io, node)


