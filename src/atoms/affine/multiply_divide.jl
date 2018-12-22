#############################################################################
# multiply_divide.jl
# Handles scalar multiplication, matrix multiplication, and scalar division
# of variables, constants and expressions.
# All expressions and atoms are subtpyes of AbstractExpr.
# Please read expressions.jl first.
#############################################################################

import Base.Broadcast.broadcasted
export sign, monotonicity, curvature, evaluate, conic_form!

### Scalar and matrix multiplication

struct MultiplyAtom <: AbstractExpr
    head::Symbol
    id_hash::UInt64
    children::Tuple{AbstractExpr, AbstractExpr}
    size::Tuple{Int, Int}

    function MultiplyAtom(x::AbstractExpr, y::AbstractExpr)
        if x.size == (1, 1)
            sz = y.size
        elseif y.size == (1, 1)
            sz = x.size
        elseif x.size[2] ==  y.size[1]
            sz = (x.size[1], y.size[2])
        else
            error("Cannot multiply two expressions of sizes $(x.size) and $(y.size)")
        end
        children = (x, y)
        return new(:*, hash(children), children, sz)
    end
end

function sign(x::MultiplyAtom)
    return sign(x.children[1]) * sign(x.children[2])
end

function monotonicity(x::MultiplyAtom)
    return (sign(x.children[2]) * Nondecreasing(), sign(x.children[1]) * Nondecreasing())
end

# Multiplication has an indefinite hessian, so if neither children are constants,
# the curvature of the atom will violate DCP.
function curvature(x::MultiplyAtom)
    if vexity(x.children[1]) != ConstVexity() && vexity(x.children[2]) != ConstVexity()
        return NotDcp()
    else
        return ConstVexity()
    end
end

function evaluate(x::MultiplyAtom)
    return evaluate(x.children[1]) * evaluate(x.children[2])
end

function conic_form!(x::MultiplyAtom, unique_conic_forms::UniqueConicForms=UniqueConicForms())
    if !has_conic_form(unique_conic_forms, x)
        # scalar multiplication
        if x.children[1].size == (1, 1) || x.children[2].size == (1, 1)
            if vexity(x.children[1]) == ConstVexity()
                const_child = x.children[1]
                expr_child = x.children[2]
            elseif vexity(x.children[2]) == ConstVexity()
                const_child = x.children[2]
                expr_child = x.children[1]
            else
                error("multiplication of two non-constant expressions is not DCP compliant")
            end
            objective = conic_form!(expr_child, unique_conic_forms)

            # make sure all 1x1 sized objects are interpreted as scalars, since
            # [1] * [1, 2, 3] is illegal in julia, but 1 * [1, 2, 3] is ok
            if const_child.size == (1, 1)
                const_multiplier = evaluate(const_child)[1]
            else
                const_multiplier = reshape(evaluate(const_child), length(const_child), 1)
            end

            objective = const_multiplier * objective

        # left matrix multiplication
        elseif x.children[1].head == :constant
            objective = conic_form!(x.children[2], unique_conic_forms)
            objective = Kron(Eye{Float64}(size(x, 2)), x.children[1].value) * objective
        # right matrix multiplication
        else
            objective = conic_form!(x.children[1], unique_conic_forms)
            objective = Kron(x.children[2].value', Eye{Float64}(size(x, 1))) * objective
        end
        cache_conic_form!(unique_conic_forms, x, objective)
    end
    return get_conic_form(unique_conic_forms, x)
end

function *(x::AbstractExpr, y::AbstractExpr)
    if hash(x) == hash(y) && x.size == (1, 1)
        return square(x)
    end
    return MultiplyAtom(x, y)
end

*(x::Value, y::AbstractExpr) = MultiplyAtom(Constant(x), y)
*(x::AbstractExpr, y::Value) = MultiplyAtom(x, Constant(y))
/(x::AbstractExpr, y::Value) = MultiplyAtom(x, Constant(1 ./ y))

### .*
# All constructors of this check (and so this function requires)
# that the first child be constant to have the expression be DCP
struct DotMultiplyAtom <: AbstractExpr
    head::Symbol
    id_hash::UInt64
    children::Tuple{AbstractExpr, AbstractExpr}
    size::Tuple{Int, Int}

    function DotMultiplyAtom(x::AbstractExpr, y::AbstractExpr)
        # check that the sizes of x and y are compatible
        try
            ones(size(x)) .* ones(size(y))
        catch
            error("cannot compute $x .* $y: sizes are not compatible")
        end
        children = (x, y)
        return new(:.*, hash(children), children, y.size)
    end
end

function sign(x::DotMultiplyAtom)
    return sign(x.children[1]) * sign(x.children[2])
end

function monotonicity(x::DotMultiplyAtom)
    return (sign(x.children[2]) * Nondecreasing(), sign(x.children[1]) * Nondecreasing())
end

function curvature(x::DotMultiplyAtom)
    if vexity(x.children[1]) == ConstVexity()
        return ConstVexity()
    else
        return NotDcp()
    end
end

function evaluate(x::DotMultiplyAtom)
    return evaluate(x.children[1]) .* evaluate(x.children[2])
end

function conic_form!(x::DotMultiplyAtom, unique_conic_forms::UniqueConicForms=UniqueConicForms())
    if !has_conic_form(unique_conic_forms, x)
        if vexity(x.children[1]) != ConstVexity()
            if vexity(x.children[2]) != ConstVexity()
                error("multiplication of two non-constant expressions is not DCP compliant")
            else
                # make sure first child is the one that's constant
                x.children[1], x.children[2] = x.children[2], x.children[1]
            end
        end
        # promote the size of the coefficient matrix, so eg
        # 3 .* x
        # works regardless of the size of x
        coeff = x.children[1].value .* Ones(size(x.children[2]))
        # promote the size of the variable
        # we've previously ensured neither x nor y is 1x1
        # and that the sizes are compatible,
        # so if the sizes aren't equal the smaller one is size 1
        var = x.children[2]
        if size(var, 1) < size(coeff, 1)
            var = Ones(size(coeff, 1)) * var
        elseif size(var, 2) < size(coeff, 2)
            var = var * Ones(1, size(coeff, 1))
        end

        const_multiplier = spdiagm(0 => vec(coeff))
        objective = const_multiplier * conic_form!(var, unique_conic_forms)
        cache_conic_form!(unique_conic_forms, x, objective)
    end
    return get_conic_form(unique_conic_forms, x)
end

function broadcasted(::typeof(*), x::Constant, y::AbstractExpr)
    if x.size == (1, 1) || y.size == (1, 1)
        return x * y
    elseif size(y, 1) < size(x, 1) && size(y, 1) == 1
        return DotMultiplyAtom(x, Ones(size(x, 1)) * y)
    elseif size(y, 2) < size(x, 2) && size(y, 2) == 1
        return DotMultiplyAtom(x, y * Ones(1, size(x, 1)))
    else
        return DotMultiplyAtom(x, y)
    end
end
broadcasted(::typeof(*), y::AbstractExpr, x::Constant) = DotMultiplyAtom(x, y)

# if neither is a constant it's not DCP, but might be nice to support anyway for eg MultiConvex
function broadcasted(::typeof(*), x::AbstractExpr, y::AbstractExpr)
    if x.size == (1, 1) || y.size == (1, 1)
        return x * y
    elseif vexity(x) == ConstVexity()
        return DotMultiplyAtom(x, y)
    elseif hash(x) == hash(y)
        return square(x)
    else
        return DotMultiplyAtom(y, x)
    end
end
broadcasted(::typeof(*), x::Value, y::AbstractExpr) = DotMultiplyAtom(Constant(x), y)
broadcasted(::typeof(*), x::AbstractExpr, y::Value) = DotMultiplyAtom(Constant(y), x)
broadcasted(::typeof(/), x::AbstractExpr, y::Value) = DotMultiplyAtom(Constant(1 ./ y), x)
# x ./ y and x / y for x constant, y variable is defined in second_order_cone.qol_elemwise.jl
