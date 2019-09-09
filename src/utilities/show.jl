import Base.show, Base.summary
export show, summary
using AbstractTrees: AbstractTrees
using .TreePrint

"""
    const MAXDEPTH = Ref(3)

Controls depth of tree printing globally for Convex.jl
"""
const MAXDEPTH = Ref(3)

"""
    show_id(io::IO, x::Union{AbstractExpr, Constraint}; digits = 4)

Print a truncated version of the objects `id_hash` field.

## Example

```julia-repl
julia> x = Variable();

julia> Convex.show_id(stdout, x)
id: 6201...
```
"""
show_id(io::IO, x::Union{AbstractExpr, Constraint}; digits = 4) = print(io, show_id(x))

function show_id(x::Union{AbstractExpr, Constraint}; digits = 4)
    hash_str = string(x.id_hash)
    hash_str = hash_str[1:nextind(hash_str, digits-1)]
    return "id: " * hash_str * "..."
end

"""
    Base.summary(io::IO, x::Variable)

Prints a one-line summary of a variable `x` to `io`.

## Examples
```julia-repl
julia> x = ComplexVariable(3,2);

julia> summary(stdout, x)
3×2 complex variable (id: 5455...)
```
"""
function Base.summary(io::IO, x::Variable)
    sgn = summary(sign(x))
    cst = vexity(x) == ConstVexity() ? " (fixed)" : ""
    cst = cst * " (" * sprint(show_id, x) * ")"
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
# julia> Variable(3,4)
# Variable
# size: (3, 4)
# sign: real
# vexity: affine
# id: 7385...
# here, the `id` will change from run to run.
function show(io::IO, x::Variable)
    print(io, "Variable")
    print(io, "\nsize: $(size(x))")
    print(io, "\nsign: ")
    summary(io, sign(x))
    print(io, "\nvexity: ")
    summary(io, vexity(x))
    println(io)
    show_id(io, x)
    if x.value !== nothing
        print(io, "\nvalue: $(x.value)")
    end
end

"""
    print_tree_rstrip(io::IO, x)

Prints the results of `TreePrint.print_tree(io, x)`
without the final newline. Used for `show` methods which
invoke `print_tree`.
"""
function print_tree_rstrip(io::IO, x)
    str = sprint(TreePrint.print_tree, x, MAXDEPTH[])
    print(io, rstrip(str))
end

# This object is used to work around the fact that
# Convex overloads booleans for AbstractExpr's
# in order to generate constraints. This is problematic
# for `AbstractTrees.print_tree` which wants to compare
# the root of the tree to itself at some point.
# By wrapping all tree roots in structs, this comparison
# occurs on the level of the `struct`, and `==` falls
# back to object equality (`===`), which is what we
# want in this case.
#
# The same construct is used below for other tree roots.
struct ConstraintRoot
    constraint::Constraint
end

TreePrint.print_tree(io::IO, c::Constraint, maxdepth = 5) = TreePrint.print_tree(io, ConstraintRoot(c), maxdepth)
AbstractTrees.children(c::ConstraintRoot) = AbstractTrees.children(c.constraint)
AbstractTrees.printnode(io::IO, c::ConstraintRoot) = AbstractTrees.printnode(io, c.constraint)

show(io::IO, c::Constraint) = print_tree_rstrip(io, c)

struct ExprRoot
    expr::AbstractExpr
end
TreePrint.print_tree(io::IO, e::AbstractExpr, maxdepth = 5) = TreePrint.print_tree(io, ExprRoot(e), maxdepth)
AbstractTrees.children(e::ExprRoot) = AbstractTrees.children(e.expr)
AbstractTrees.printnode(io::IO, e::ExprRoot) = AbstractTrees.printnode(io, e.expr)


show(io::IO, e::AbstractExpr) = print_tree_rstrip(io, e)


struct ProblemObjectiveRoot
    head::Symbol
    objective::AbstractExpr
end

AbstractTrees.children(p::ProblemObjectiveRoot) = (p.objective,)
AbstractTrees.printnode(io::IO, p::ProblemObjectiveRoot) =  print(io, string(p.head))

struct ProblemConstraintsRoot
    constraints::Vector{Constraint}
end

AbstractTrees.children(p::ProblemConstraintsRoot) = p.constraints
AbstractTrees.printnode(io::IO, p::ProblemConstraintsRoot) = print(io, "subject to")


function TreePrint.print_tree(io::IO, p::Problem, maxdepth = 5)
    TreePrint.print_tree(io, ProblemObjectiveRoot(p.head, p.objective), maxdepth)
    if !(isempty(p.constraints))
        TreePrint.print_tree(io, ProblemConstraintsRoot(p.constraints), maxdepth)
    end
end

function show(io::IO, p::Problem)
    TreePrint.print_tree(io, p, MAXDEPTH[])
    print(io, "\ncurrent status: $(p.status)")
    if p.status == "solved"
        print(io, " with optimal value of $(round(p.optval, digits=4))")
    end
end
