# This module is needed until AbstractTrees.jl#37 is fixed.
# (PR: https://github.com/Keno/AbstractTrees.jl/pull/38)
# because currently `print_tree` does not respect `maxdepth`.
# This just implements the changes in the above PR.
# Code in this file is modified from AbstractTrees.jl
# See LICENSE for a copy of its MIT license.
module TreePrint

using AbstractTrees: printnode, treekind, IndexedTree, children

# Printing
struct TreeCharSet
    mid::Any
    terminator::Any
    skip::Any
    dash::Any
    ellipsis::Any
end

# Default charset
TreeCharSet() = TreeCharSet('├', '└', '│', '─', '…')

function print_prefix(io, depth, charset, active_levels)
    for current_depth in 0:(depth-1)
        if current_depth in active_levels
            print(io, charset.skip, " "^(textwidth(charset.dash) + 1))
        else
            print(
                io,
                " "^(textwidth(charset.skip) + textwidth(charset.dash) + 1),
            )
        end
    end
end

@doc raw"""
# Usage
Prints an ASCII formatted representation of the `tree` to the given `io` object.
By default all children will be printed up to a maximum level of 5, though this
valud can be overriden by the `maxdepth` parameter. The charset to use in
printing can be customized using the `charset` keyword argument.

# Examples
```julia
julia> print_tree(STDOUT,Dict("a"=>"b","b"=>['c','d']))
Dict{String,Any}("b"=>['c','d'],"a"=>"b")
├─ b
│  ├─ c
│  └─ d
└─ a
   └─ b

julia> print_tree(STDOUT,Dict("a"=>"b","b"=>['c','d']);
        charset = TreeCharSet('+','\\','|',"--"))
Dict{String,Any}("b"=>['c','d'],"a"=>"b")
+-- b
|   +-- c
|   \-- d
\-- a
   \-- b
```

"""
print_tree

function _print_tree(
    printnode::Function,
    io::IO,
    tree,
    maxdepth = 5,
    maxwidth = Inf;
    depth = 0,
    active_levels = Int[],
    charset = TreeCharSet(),
    withinds = false,
    inds = [],
    from = nothing,
    to = nothing,
    roottree = tree,
)
    nodebuf = IOBuffer()
    isa(io, IOContext) && (nodebuf = IOContext(nodebuf, io))
    if withinds
        printnode(nodebuf, tree, inds)
    else
        tree != roottree && isa(treekind(roottree), IndexedTree) ?
        printnode(nodebuf, roottree[tree]) : printnode(nodebuf, tree)
    end
    str = String(take!(isa(nodebuf, IOContext) ? nodebuf.io : nodebuf))
    for (i, line) in enumerate(split(str, '\n'))
        i != 1 && print_prefix(io, depth, charset, active_levels)
        println(io, line)
    end
    depth > maxdepth && return
    c =
        isa(treekind(roottree), IndexedTree) ? childindices(roottree, tree) :
        children(roottree, tree)
    if c !== ()
        width = 0
        s = Iterators.Stateful(
            from === nothing ? pairs(c) : Iterators.Rest(pairs(c), from),
        )
        while !isempty(s) && width < maxwidth
            width += 1
            ind, child = popfirst!(s)
            ind === to && break
            active = false
            child_active_levels = active_levels
            print_prefix(io, depth, charset, active_levels)
            if isempty(s)
                print(io, charset.terminator)
            else
                print(io, charset.mid)
                child_active_levels = push!(copy(active_levels), depth)
            end
            print(io, charset.dash, ' ')
            print_tree(
                depth == maxdepth ? (io, val) -> print(io, charset.ellipsis) :
                printnode,
                io,
                child,
                maxdepth,
                maxwidth;
                depth = depth + 1,
                active_levels = child_active_levels,
                charset = charset,
                withinds = withinds,
                inds = withinds ? [inds; ind] : [],
                roottree = roottree,
            )
        end
        if !isempty(s)
            print_prefix(io, depth, charset, active_levels)
            println(io, "⋮")
        end
    end
end
function print_tree(f::Function, io::IO, tree, args...; kwargs...)
    return _print_tree(f, io, tree, args...; kwargs...)
end
function print_tree(io::IO, tree, args...; kwargs...)
    return print_tree(printnode, io, tree, args...; kwargs...)
end
function print_tree(tree, args...; kwargs...)
    return print_tree(stdout::IO, tree, args...; kwargs...)
end
end
