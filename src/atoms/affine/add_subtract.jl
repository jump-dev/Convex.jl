#############################################################################
# add_subtract.jl
# Handles unary negation, addition and subtraction of variables, constants
# and expressions.
# All expressions and atoms are subtpyes of AbstractExpr.
# Please read expressions.jl first.
#############################################################################

export +, -
export sign, curvature, monotonicity, evaluate

### Unary Negation

struct NegateAtom <: AbstractExpr
    head::Symbol
    id_hash::UInt64
    children::Tuple{AbstractExpr}
    size::Tuple{Int, Int}

    function NegateAtom(x::AbstractExpr)
        children = (x,)
        return new(:-, hash(children), children, x.size)
    end
end

function sign(x::NegateAtom)
    return -sign(x.children[1])
end

function monotonicity(x::NegateAtom)
    return (Nonincreasing(),)
end

function curvature(x::NegateAtom)
    return ConstVexity()
end

function evaluate(x::NegateAtom)
    return -evaluate(x.children[1])
end

-(x::AbstractExpr) = NegateAtom(x)

function conic_form!(x::NegateAtom, unique_conic_forms::UniqueConicForms=UniqueConicForms())
    if !has_conic_form(unique_conic_forms, x)
        objective = conic_form!(x.children[1], unique_conic_forms)
        objective = -objective
        cache_conic_form!(unique_conic_forms, x, objective)
    end
    return get_conic_form(unique_conic_forms, x)
end


### Addition
struct AdditionAtom <: AbstractExpr
    head::Symbol
    id_hash::UInt64
    children::Array{AbstractExpr, 1}
    size::Tuple{Int, Int}

    function AdditionAtom(x::AbstractExpr, y::AbstractExpr)
        # find the size of the expression = max of size of x and size of y
        if x.size == y.size || y.size == (1, 1)
            sz = x.size
        elseif x.size == (1, 1)
            sz = y.size
        else
            error("cannot add expressions of sizes $(x.size) and $(y.size)")
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
        return new(:+, hash(children), children, sz)
    end
end

function sign(x::AdditionAtom)
    return sum(Sign[sign(child) for child in x.children])
    # Creating an array of type Sign and adding all the sign of xhildren of x so if anyone is complex the resultant sign would be complex.
end

function monotonicity(x::AdditionAtom)
    return Monotonicity[Nondecreasing() for child in x.children]
end

function curvature(x::AdditionAtom)
    return ConstVexity()
end

function evaluate(x::AdditionAtom)
    # broadcast is used in reduction instead of using sum directly to support addition
    # between scalars and arrays
    return mapreduce(evaluate, (a, b) -> a .+ b, x.children)
end

function conic_form!(x::AdditionAtom, unique_conic_forms::UniqueConicForms=UniqueConicForms())
    if !has_conic_form(unique_conic_forms, x)
        objective = ConicObj()
        for child in x.children
            child_objective = conic_form!(child, unique_conic_forms)
            if x.size != child.size
                child_objective = promote_size(child_objective, length(x))
            end
            objective += child_objective
        end
        cache_conic_form!(unique_conic_forms, x, objective)
    end
    return get_conic_form(unique_conic_forms, x)
end

+(x::AbstractExpr, y::AbstractExpr) = AdditionAtom(x, y)
+(x::Value, y::AbstractExpr) = AdditionAtom(Constant(x), y)
+(x::AbstractExpr, y::Value) = AdditionAtom(x, Constant(y))
-(x::AbstractExpr, y::AbstractExpr) = x + (-y)
-(x::Value, y::AbstractExpr) = Constant(x) + (-y)
-(x::AbstractExpr, y::Value) = x + Constant(-y)
