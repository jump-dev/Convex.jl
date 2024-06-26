# The `TreePrint` module contains code from the AbstractTrees.jl package
# (https://github.com/Keno/AbstractTrees.jl) which is licensed under the
# MIT "Expat" License:
#
# This module originally existed for AbstractTrees.jl#37, but has since diverged
# from the functionality of AbstractTrees.print_tree. It is now a separate
# implementation of tree printing.
#
# > Copyright (c) 2015: Keno Fischer.
# >
# > Permission is hereby granted, free of charge, to any person obtaining
# > a copy of this software and associated documentation files (the
# > "Software"), to deal in the Software without restriction, including
# > without limitation the rights to use, copy, modify, merge, publish,
# > distribute, sublicense, and/or sell copies of the Software, and to
# > permit persons to whom the Software is furnished to do so, subject to
# > the following conditions:
# >
# > The above copyright notice and this permission notice shall be
# > included in all copies or substantial portions of the Software.
# >
# > THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# > EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# > MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
# > IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
# > CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
# > TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
# > SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

module TreePrint

import AbstractTrees

struct TreeCharSet
    mid::Any
    terminator::Any
    skip::Any
    dash::Any
    ellipsis::Any
end

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
    return
end

"""
    print_tree(
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
|   \\-- d
\\-- a
   \\-- b
```
"""
function print_tree(
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
    if io isa IOContext
        nodebuf = IOContext(nodebuf, io)
    end
    if withinds
        printnode(nodebuf, tree, inds)
    else
        printnode(nodebuf, tree)
    end
    str = String(take!(isa(nodebuf, IOContext) ? nodebuf.io : nodebuf))
    for (i, line) in enumerate(split(str, '\n'))
        if i != 1
            print_prefix(io, depth, charset, active_levels)
        end
        println(io, line)
    end
    if depth > maxdepth
        return
    end
    c = AbstractTrees.children(tree)
    if isempty(c)
        return
    end
    width = 0
    s = Iterators.Stateful(
        from === nothing ? pairs(c) : Iterators.Rest(pairs(c), from),
    )
    while !isempty(s) && width < maxwidth
        width += 1
        ind, child = popfirst!(s)
        if ind === to
            break
        end
        child_active_levels = active_levels
        print_prefix(io, depth, charset, active_levels)
        if isempty(s)
            print(io, charset.terminator)
        else
            print(io, charset.mid)
            child_active_levels = vcat(active_levels, depth)
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
    return
end

function print_tree(io::IO, tree, args...; kwargs...)
    return print_tree(AbstractTrees.printnode, io, tree, args...; kwargs...)
end

function print_tree(tree, args...; kwargs...)
    return print_tree(stdout, tree, args...; kwargs...)
end

end
