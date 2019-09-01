import Base.show, Base.summary
export show, summary
import AbstractTrees


function Base.summary(io::IO, x::Variable)
    hash_str = string(x.id_hash)
    hash_str = hash_str[1:nextind(hash_str, 3)]
    sgn = summary(sign(x))
    cst = vexity(x) == ConstVexity() ? " (fixed)" : ""
    cst = cst * " (id: " * hash_str * "...)"
    if size(x) == (1,1)
        print(io, "$(sgn) variable$(cst)")
    elseif size(x,2) == 1
        print(io, "$(size(x,1))-element $(sgn) variable$(cst)")
    else
        print(io, "$(size(x,1))×$(size(x,2)) $(sgn) variable$(cst)")

    end
end

Base.summary(io::IO, ::AffineVexity) = print(io, "affine")
Base.summary(io::IO, ::ConvexVexity) = print(io, "convex")
Base.summary(io::IO, ::ConcaveVexity) = print(io, "concave")
Base.summary(io::IO, ::ConstVexity) = print(io, "constant")

Base.summary(io::IO, ::Positive) = print(io, "positive")
Base.summary(io::IO, ::Negative) = print(io, "negative")
Base.summary(io::IO, ::NoSign) = print(io, "real")
Base.summary(io::IO, ::ComplexSign) = print(io, "complex")

function Base.summary(io::IO, c::Constraint)
    print(io, "$(c.head) constraint (")
    summary(io, vexity(c))
    print(io, ")")
end

function Base.summary(io::IO, e::AbstractExpr)
    print(io, "$(e.head) (")
    summary(io, vexity(e))
    print(io, "; ")
    summary(io, sign(e))
    print(io, ")")
end

# A Constant is simply a wrapper around a native Julia constant
# Hence, we simply display its value
show(io::IO, x::Constant) = print(io, x.value)

# A variable, for example, Variable(3, 4), will be displayed as:
# Variable of
# size: (3, 4)
# sign: NoSign()
# vexity: AffineVexity()
function show(io::IO, x::Variable)
    print(io, """Variable of
          size: ($(x.size[1]), $(x.size[2]))
          sign: $(x.sign)
          vexity: $(x.vexity)""")
    if x.value !== nothing
        print(io, "\nvalue: $(x.value)")
    end
end

struct ShowConstraint
    constraint::Constraint
end

AbstractTrees.children(c::ShowConstraint) = AbstractTrees.children(c.constraint)
AbstractTrees.printnode(io::IO, c::ShowConstraint) = AbstractTrees.printnode(io, c.constraint)

show(io::IO, c::Constraint) = AbstractTrees.print_tree(io, ShowConstraint(c))



struct ShowExpr
    expr::AbstractExpr
end

AbstractTrees.children(c::ShowExpr) = AbstractTrees.children(c.expr)
AbstractTrees.printnode(io::IO, c::ShowExpr) = AbstractTrees.printnode(io, c.expr)

show(io::IO, c::AbstractExpr) = AbstractTrees.print_tree(io, ShowExpr(c))


struct ShowProblemObjective
    head::Symbol
    objective::AbstractExpr
end

AbstractTrees.children(p::ShowProblemObjective) = (p.objective,)
AbstractTrees.printnode(io::IO, p::ShowProblemObjective) =  print(io,string(p.head))

struct ShowProblemConstraints
    constraints::Vector{Constraint}
end

AbstractTrees.children(p::ShowProblemConstraints) = tuple(p.constraints)
AbstractTrees.printnode(io::IO, p::ShowProblemConstraints) = print(io,"subject to")


function show(io::IO, p::Problem)
    AbstractTrees.print_tree(io, ShowProblemObjective(p.head, p.objective))
    AbstractTrees.print_tree(io, ShowProblemConstraints(p.constraints))
    print(io, "\ncurrent status: $(p.status)")
    if p.status == "solved"
        print(io, " with optimal value of $(round(p.optval, digits=4))")
    end
end
