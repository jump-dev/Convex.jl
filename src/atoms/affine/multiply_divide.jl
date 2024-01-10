mutable struct MultiplyAtom <: AbstractExpr
    children::Tuple{AbstractExpr,AbstractExpr}
    size::Tuple{Int,Int}

    function MultiplyAtom(x::AbstractExpr, y::AbstractExpr)
        sz = if x.size == (1, 1)
            y.size
        elseif y.size == (1, 1)
            x.size
        elseif x.size[2] == y.size[1]
            (x.size[1], y.size[2])
        else
            error(
                "Cannot multiply two expressions of sizes $(x.size) and $(y.size)",
            )
        end
        return new((x, y), sz)
    end
end

head(io::IO, ::MultiplyAtom) = print(io, "*")

Base.sign(x::MultiplyAtom) = sign(x.children[1]) * sign(x.children[2])

function monotonicity(x::MultiplyAtom)
    return (
        sign(x.children[2]) * Nondecreasing(),
        sign(x.children[1]) * Nondecreasing(),
    )
end

# Multiplication has an indefinite hessian, so if neither children are constants,
# the curvature of the atom will violate DCP.
function curvature(x::MultiplyAtom)
    if vexity(x.children[1]) != ConstVexity() &&
       vexity(x.children[2]) != ConstVexity()
        return NotDcp()
    end
    return ConstVexity()
end

evaluate(x::MultiplyAtom) = evaluate(x.children[1]) * evaluate(x.children[2])

complex_convert(::Type{T}, x) where {T} = real_convert(Complex{T}, x)

real_convert(::Type{T}, x::Number) where {T} = T(x)

real_convert(::Type{T}, x::AbstractMatrix) where {T} = create_sparse(T, x)

real_convert(::Type{T}, x::SPARSE_MATRIX{T}) where {T} = x

real_convert(::Type{T}, x::SPARSE_VECTOR{T}) where {T} = x

real_convert(::Type{T}, x::AbstractVector) where {T} = SPARSE_VECTOR{T}(x)

function new_conic_form!(context::Context{T}, x::MultiplyAtom) where {T}
    if x.children[1].size == (1, 1) || x.children[2].size == (1, 1)
        # scalar multiplication
        if vexity(x.children[1]) == ConstVexity()
            const_child, expr_child = x.children
        elseif vexity(x.children[2]) == ConstVexity()
            expr_child, const_child = x.children
        else
            error(
                "multiplication of two non-constant expressions is not DCP compliant",
            )
        end
        objective = conic_form!(context, expr_child)
        # make sure all 1x1 sized objects are interpreted as scalars, since
        # [1] * [1, 2, 3] is illegal in julia, but 1 * [1, 2, 3] is ok
        const_multiplier = if const_child.size == (1, 1)  # Scalar
            evaluate(const_child)[1]
        else  # Matrix
            reshape(evaluate(const_child), length(const_child), 1)
        end
        const_multiplier = if iscomplex(const_multiplier)
            complex_convert(T, const_multiplier)
        else
            real_convert(T, const_multiplier)
        end
        return operate(add_operation, T, sign(x), const_multiplier, objective)
    elseif vexity(x.children[1]) == ConstVexity()
        # left matrix multiplication
        objective = conic_form!(context, x.children[2])
        const_multiplier = evaluate(x.children[1])
        const_multiplier = if iscomplex(const_multiplier)
            complex_convert(T, const_multiplier)
        else
            real_convert(T, const_multiplier)
        end
        return operate(
            add_operation,
            T,
            sign(x),
            kron(spidentity(T, x.size[2]), const_multiplier),
            objective,
        )
    else
        # right matrix multiplication
        objective = conic_form!(context, x.children[1])
        const_multiplier = evaluate(x.children[2])
        const_multiplier = if iscomplex(const_multiplier)
            complex_convert(T, const_multiplier)
        else
            real_convert(T, const_multiplier)
        end
        return operate(
            add_operation,
            T,
            sign(x),
            kron(transpose(const_multiplier), spidentity(T, x.size[1])),
            objective,
        )
    end
end

function Base.:*(x::AbstractExpr, y::AbstractExpr)
    if isequal(x, y) && x.size == (1, 1)
        return square(x)
    end
    return MultiplyAtom(x, y)
end

Base.:*(x::Value, y::AbstractExpr) = MultiplyAtom(constant(x), y)

Base.:*(x::AbstractExpr, y::Value) = MultiplyAtom(x, constant(y))

Base.:/(x::AbstractExpr, y::Value) = MultiplyAtom(x, constant(1 ./ y))

# ambiguity
#
# function Base.:(*)(
#     x::Convex.AbstractExpr,
#     y::Union{
#         LinearAlgebra.Transpose{
#             <:Any,
#             <:SuiteSparseGraphBLAS.AbstractGBArray{T,F,O},
#         },
#         SuiteSparseGraphBLAS.AbstractGBArray{T,F,O},
#     },
# ) where {T,F,O}
#     return MultiplyAtom(x, constant(y))
# end
#
# function Base.:(*)(
#     x::Union{
#         LinearAlgebra.Transpose{
#             <:Any,
#             <:SuiteSparseGraphBLAS.AbstractGBArray{T,F,O},
#         },
#         SuiteSparseGraphBLAS.AbstractGBArray{T,F,O},
#     },
#     y::Convex.AbstractExpr,
# ) where {T,F,O}
#     return MultiplyAtom(constant(x), y)
# end

function dotmultiply(x, y)
    if size(x) == (1, 1) || size(y) == (1, 1)
        return x * y
    end
    if vexity(x) != ConstVexity()
        if vexity(y) != ConstVexity()
            error(
                "multiplication of two non-constant expressions is not DCP compliant",
            )
        end
        x, y = y, x
    end
    # promote the size of the coefficient matrix, so e.g., 3 .* x works
    # regardless of the size of x
    coeff = evaluate(x) .* ones(size(y))
    # promote the size of the variable
    # we've previously ensured neither x nor y is 1x1
    # and that the sizes are compatible,
    # so if the sizes aren't equal the smaller one is size 1
    var = y
    if size(var, 1) < size(coeff, 1)
        var = ones(size(coeff, 1)) * var
    elseif size(var, 2) < size(coeff, 2)
        var = var * ones(1, size(coeff, 1))
    end
    const_multiplier = LinearAlgebra.Diagonal(vec(coeff))
    return reshape(const_multiplier * vec(var), size(var)...)
end

# if neither is a constant it's not DCP, but might be nice to support anyway for
# eg MultiConvex
function Base.Broadcast.broadcasted(
    ::typeof(*),
    x::AbstractExpr,
    y::AbstractExpr,
)
    if isequal(x, y)
        return square(x)
    end
    return dotmultiply(x, y)
end

function Base.Broadcast.broadcasted(::typeof(*), x::Value, y::AbstractExpr)
    return dotmultiply(constant(x), y)
end

function Base.Broadcast.broadcasted(::typeof(*), x::AbstractExpr, y::Value)
    return dotmultiply(constant(y), x)
end

function Base.Broadcast.broadcasted(::typeof(/), x::AbstractExpr, y::Value)
    return dotmultiply(constant(1 ./ y), x)
end

# x ./ y and x / y for x constant, y variable is defined in
# second_order_cone/qol_elemwise.jl
