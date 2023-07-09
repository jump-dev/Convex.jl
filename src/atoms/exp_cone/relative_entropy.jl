#############################################################################
# relative_entropy.jl
# relative entropy (ie, sum_i( x_i log (x_i/y_i) ) of expressions x and y
# All expressions and atoms are subtypes of AbstractExpr.
# Please read expressions.jl first.
#############################################################################

# TODO: make this work for a *list* of inputs, rather than just for scalar/vector/matrix inputs

struct RelativeEntropyAtom <: AbstractExpr
    head::Symbol
    id_hash::UInt64
    children::Tuple{AbstractExpr,AbstractExpr}
    size::Tuple{Int,Int}

    function RelativeEntropyAtom(x::AbstractExpr, y::AbstractExpr)
        if sign(x) == ComplexSign() || sign(y) == ComplexSign()
            error(
                "Both the arguments should be real but these are instead $(sign(x)) and $(sign(y))",
            )
        else
            children = (x, y)
            return new(:relative_entropy, hash(children), children, (1, 1))
        end
    end
end

function sign(x::RelativeEntropyAtom)
    return NoSign()
end

function monotonicity(x::RelativeEntropyAtom)
    return (NoMonotonicity(), NoMonotonicity())
end

function curvature(x::RelativeEntropyAtom)
    return ConvexVexity()
end

function evaluate(e::RelativeEntropyAtom)
    x = vectorize(evaluate(e.children[1]))
    y = vectorize(evaluate(e.children[2]))
    if any(isnan, y)
        return Inf
    end

    out = x .* log.(x ./ y)
    # fix value when x=0:
    # out will only be NaN if x=0, in which case the correct value is 0
    out[isnan.(out)] .= 0
    return sum(out)
end

function conic_form!(context::Context{T}, e::RelativeEntropyAtom) where {T}
    # relative_entropy(x,y) = sum_i( x_i log (x_i/y_i) )
    w = conic_form!(context, e.children[1])
    v = conic_form!(context, e.children[2])
    u = conic_form!(context, Variable())
    f = operate(vcat, T, sign(e), u, v, w)
    d = MOI.output_dimension(w)
    @assert d == MOI.output_dimension(v)
    MOI_add_constraint(context.model, f, MOI.RelativeEntropyCone(2d + 1))
    return u
end

function test(optimizer)
    model = Convex.Context{Float64}(optimizer).model
    u = MOI.add_variable(model)
    v = MOI.add_variable(model)
    w = MOI.add_variable(model)
    # f = MOI.VectorOfVariables([u, v, w])
    f = MOI.VectorAffineFunction{Float64}(
        [
            MOI.VectorAffineTerm{Float64}(
                1,
                MOI.ScalarAffineTerm{Float64}(1.0, u),
            ),
            MOI.VectorAffineTerm{Float64}(
                2,
                MOI.ScalarAffineTerm{Float64}(1.0, v),
            ),
            MOI.VectorAffineTerm{Float64}(
                3,
                MOI.ScalarAffineTerm{Float64}(1.0, w),
            ),
        ],
        [0.0, 0.0, 0.0],
    )
    MOI.add_constraint(model, f, MOI.RelativeEntropyCone(1))

    g2 = MOI.VectorAffineFunction{Float64}(
        [
            MOI.VectorAffineTerm{Float64}(
                1,
                MOI.ScalarAffineTerm{Float64}(1.0, v),
            ),
        ],
        [-5.0],
    )
    MOI.add_constraint(model, g2, MOI.Zeros(1))

    g1 = MOI.VectorAffineFunction{Float64}(
        [
            MOI.VectorAffineTerm{Float64}(
                1,
                MOI.ScalarAffineTerm{Float64}(1.0, w),
            ),
        ],
        [-10.0],
    )
    MOI.add_constraint(model, g1, MOI.Nonpositives(1))
    # MOI.add_constraint(model, MOI.SingleVariable(v), MOI.EqualTo(5.0))
    # MOI.add_constraint(model, MOI.SingleVariable(w), MOI.LessThan(10.0))
    # o=MOI.SingleVariable(u)
    o = MOI.ScalarAffineFunction{Float64}(
        MOI.ScalarAffineTerm{Float64}[MOI.ScalarAffineTerm{Float64}(-1.0, u)],
        0.0,
    )
    MOI.set(model, MOI.ObjectiveFunction{typeof(o)}(), o)
    MOI.set(model, MOI.ObjectiveSense(), MOI.MAX_SENSE)
    MOI.optimize!(model)
    obj_val = MOI.get(model, MOI.ObjectiveValue())
    u_val = MOI.get(model, MOI.VariablePrimal(), u)
    return obj_val, u_val
end

relative_entropy(x::AbstractExpr, y::AbstractExpr) = RelativeEntropyAtom(x, y)
# y*log(x/y)
log_perspective(x::AbstractExpr, y::AbstractExpr) = -relative_entropy(y, x)
