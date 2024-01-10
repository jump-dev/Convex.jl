mutable struct NegateAtom <: AbstractExpr
    children::Tuple{AbstractExpr}
    size::Tuple{Int,Int}

    NegateAtom(x::AbstractExpr) = new((x,), x.size)
end

Base.sign(x::NegateAtom) = -sign(x.children[1])

monotonicity(::NegateAtom) = (Nonincreasing(),)

curvature(::NegateAtom) = ConstVexity()

evaluate(x::NegateAtom) = -evaluate(x.children[1])

Base.:-(x::AbstractExpr) = NegateAtom(x)

Base.:-(x::Union{Constant,ComplexConstant}) = constant(-evaluate(x))

function new_conic_form!(context::Context{T}, A::NegateAtom) where {T}
    subobj = conic_form!(context, only(AbstractTrees.children(A)))
    if subobj isa Value
        return -subobj
    end
    return operate(-, T, sign(A), subobj)
end

mutable struct AdditionAtom <: AbstractExpr
    children::Array{AbstractExpr,1}
    size::Tuple{Int,Int}

    function AdditionAtom(x::AbstractExpr, y::AbstractExpr)
        # find the size of the expression = max of size of x and size of y
        if x.size == y.size || y.size == (1, 1)
            sz = x.size
            if y.size == (1, 1)
                y = y * ones(sz)
            end
        elseif x.size == (1, 1)
            sz = y.size
            x = x * ones(sz)
        else
            error("Cannot add expressions of sizes $(x.size) and $(y.size)")
        end
        if x.size != y.size
            if (x isa Constant || x isa ComplexConstant) && (x.size == (1, 1))
                x = constant(fill(evaluate(x), y.size))
            elseif (y isa Constant || y isa ComplexConstant) &&
                   (y.size == (1, 1))
                y = constant(fill(evaluate(y), x.size))
            end
        end
        # see if we're forming a sum of more than two terms and condense them
        children = AbstractExpr[]
        if isa(x, AdditionAtom)
            append!(children, x.children)
        else
            push!(children, x)
        end
        if isa(y, AdditionAtom)
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

monotonicity(x::AdditionAtom) = [Nondecreasing() for _ in x.children]

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
