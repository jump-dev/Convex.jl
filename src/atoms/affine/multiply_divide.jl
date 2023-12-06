#############################################################################
# multiply_divide.jl
# Handles scalar multiplication, matrix multiplication, and scalar division
# of variables, constants and expressions.
# All expressions and atoms are subtpyes of AbstractExpr.
# Please read expressions.jl first.
#############################################################################

import Base.Broadcast.broadcasted

### Scalar and matrix multiplication

mutable struct MultiplyAtom <: AbstractExpr
    children::Tuple{AbstractExpr,AbstractExpr}
    size::Tuple{Int,Int}

    function MultiplyAtom(x::AbstractExpr, y::AbstractExpr)
        if x.size == (1, 1)
            sz = y.size
        elseif y.size == (1, 1)
            sz = x.size
        elseif x.size[2] == y.size[1]
            sz = (x.size[1], y.size[2])
        else
            error(
                "Cannot multiply two expressions of sizes $(x.size) and $(y.size)",
            )
        end
        children = (x, y)
        return new(children, sz)
    end
end

head(io::IO, ::MultiplyAtom) = print(io, "*")

function sign(x::MultiplyAtom)
    return sign(x.children[1]) * sign(x.children[2])
end

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
    else
        return ConstVexity()
    end
end

function evaluate(x::MultiplyAtom)
    return evaluate(x.children[1]) * evaluate(x.children[2])
end

function complex_convert(::Type{T}, x) where {T}
    return real_convert(Complex{T}, x)
end
function real_convert(::Type{T}, x::Number) where {T}
    return T(x)
end
function real_convert(::Type{T}, x::AbstractMatrix) where {T}
    return create_sparse(T, x)
end

function real_convert(::Type{T}, x::SPARSE_MATRIX{T}) where {T}
    return x
end

function real_convert(::Type{T}, x::SPARSE_VECTOR{T}) where {T}
    return x
end
function real_convert(::Type{T}, x::AbstractVector) where {T}
    return SPARSE_VECTOR{T}(x)
end

function _conic_form!(context::Context{T}, x::MultiplyAtom) where {T}
    # scalar multiplication
    if x.children[1].size == (1, 1) || x.children[2].size == (1, 1)
        if vexity(x.children[1]) == ConstVexity()
            const_child = x.children[1]
            expr_child = x.children[2]
        elseif vexity(x.children[2]) == ConstVexity()
            const_child = x.children[2]
            expr_child = x.children[1]
        else
            error(
                "multiplication of two non-constant expressions is not DCP compliant",
            )
        end

        objective = conic_form!(context, expr_child)

        # make sure all 1x1 sized objects are interpreted as scalars, since
        # [1] * [1, 2, 3] is illegal in julia, but 1 * [1, 2, 3] is ok
        if const_child.size == (1, 1)
            # Scalar
            const_multiplier = evaluate(const_child)[1]
        else
            # Matrix
            const_multiplier =
                reshape(evaluate(const_child), length(const_child), 1)
        end

        if iscomplex(const_multiplier)
            const_multiplier = complex_convert(T, const_multiplier)
        else
            const_multiplier = real_convert(T, const_multiplier)
        end

        return operate(add_operation, T, sign(x), const_multiplier, objective)

        # left matrix multiplication
    elseif vexity(x.children[1]) == ConstVexity()
        objective = conic_form!(context, x.children[2])

        const_multiplier = evaluate(x.children[1])

        if iscomplex(const_multiplier)
            const_multiplier = complex_convert(T, const_multiplier)
        else
            const_multiplier = real_convert(T, const_multiplier)
        end

        return operate(
            add_operation,
            T,
            sign(x),
            kron(spidentity(T, x.size[2]), const_multiplier),
            objective,
        )

        # right matrix multiplication
    else
        objective = conic_form!(context, x.children[1])
        const_multiplier = evaluate(x.children[2])

        if iscomplex(const_multiplier)
            const_multiplier = complex_convert(T, const_multiplier)
        else
            const_multiplier = real_convert(T, const_multiplier)
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

function *(x::AbstractExpr, y::AbstractExpr)
    if isequal(x, y) && x.size == (1, 1)
        return square(x)
    end
    return MultiplyAtom(x, y)
end

*(x::Value, y::AbstractExpr) = MultiplyAtom(constant(x), y)
*(x::AbstractExpr, y::Value) = MultiplyAtom(x, constant(y))
/(x::AbstractExpr, y::Value) = MultiplyAtom(x, constant(1 ./ y))

# # ambiguity
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
    if vexity(x) != ConstVexity()
        if vexity(y) != ConstVexity()
            error(
                "multiplication of two non-constant expressions is not DCP compliant",
            )
        else
            x, y = y, x
        end
    end

    # promote the size of the coefficient matrix, so eg
    # 3 .* x
    # works regardless of the size of x
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
    const_multiplier = Diagonal(vec(coeff))
    return reshape(const_multiplier * vec(var), size(var)...)
end

function broadcasted(
    ::typeof(*),
    x::Union{Constant,ComplexConstant},
    y::AbstractExpr,
)
    if x.size == (1, 1) || y.size == (1, 1)
        return x * y
    elseif size(y, 1) < size(x, 1) && size(y, 1) == 1
        return dotmultiply(x, ones(size(x, 1)) * y)
    elseif size(y, 2) < size(x, 2) && size(y, 2) == 1
        return dotmultiply(x, y * ones(1, size(x, 1)))
    else
        return dotmultiply(x, y)
    end
end
function broadcasted(
    ::typeof(*),
    y::AbstractExpr,
    x::Union{Constant,ComplexConstant},
)
    return dotmultiply(x, y)
end

# if neither is a constant it's not DCP, but might be nice to support anyway for eg MultiConvex
function broadcasted(::typeof(*), x::AbstractExpr, y::AbstractExpr)
    if x.size == (1, 1) || y.size == (1, 1)
        return x * y
    elseif vexity(x) == ConstVexity()
        return dotmultiply(x, y)
    elseif isequal(x, y)
        return square(x)
    else
        return dotmultiply(y, x)
    end
end
function broadcasted(::typeof(*), x::Value, y::AbstractExpr)
    return dotmultiply(constant(x), y)
end
function broadcasted(::typeof(*), x::AbstractExpr, y::Value)
    return dotmultiply(constant(y), x)
end
function broadcasted(::typeof(/), x::AbstractExpr, y::Value)
    return dotmultiply(constant(1 ./ y), x)
end
# x ./ y and x / y for x constant, y variable is defined in second_order_cone.qol_elemwise.jl
