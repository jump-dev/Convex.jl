mutable struct AdditionAtom <: AbstractExpr
    children::Vector{AbstractExpr}
    size::Tuple{Int,Int}

    function AdditionAtom(x::AbstractExpr, y::AbstractExpr)
        # find the size of the expression = max of size of x and size of y
        sz = if x.size == y.size
            x.size
        elseif y.size == (1, 1)
            y = y * ones(x.size)
            x.size
        elseif x.size == (1, 1)
            x = x * ones(y.size)
            y.size
        else
            error(
                "[AdditionAtom] cannot add expressions of sizes $(x.size) and $(y.size)",
            )
        end
        # See if we're forming a sum of more than two terms and condense them
        children = AbstractExpr[]
        if x isa AdditionAtom
            append!(children, x.children)
        else
            push!(children, x)
        end
        if y isa AdditionAtom
            append!(children, y.children)
        else
            push!(children, y)
        end
        return new(children, sz)
    end
end

head(io::IO, ::AdditionAtom) = print(io, "+")

# Creating an array of type Sign and adding all the sign of children of x,
# so if anyone is complex the resultant sign would be complex.
Base.sign(x::AdditionAtom) = sum(sign.(x.children))

monotonicity(x::AdditionAtom) = ntuple(i -> Nondecreasing(), length(x.children))

curvature(::AdditionAtom) = ConstVexity()

function evaluate(x::AdditionAtom)
    # broadcast is used in reduction instead of using sum directly to support
    # addition between scalars and arrays
    return mapreduce(evaluate, (a, b) -> a .+ b, x.children)
end

function new_conic_form!(context::Context{T}, x::AdditionAtom) where {T}
    return operate(
        +,
        T,
        sign(x),
        (conic_form!(context, c) for c in AbstractTrees.children(x))...,
    )
end

Base.:+(x::AbstractExpr, y::AbstractExpr) = AdditionAtom(x, y)

Base.:+(x::Value, y::AbstractExpr) = AdditionAtom(constant(x), y)

Base.:+(x::AbstractExpr, y::Value) = AdditionAtom(x, constant(y))

Base.:-(x::AbstractExpr, y::AbstractExpr) = x + (-y)

Base.:-(x::Value, y::AbstractExpr) = constant(x) + (-y)

Base.:-(x::AbstractExpr, y::Value) = x + constant(-y)
