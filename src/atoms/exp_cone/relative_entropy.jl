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
    size::Tuple{Int, Int}

    function RelativeEntropyAtom(x::AbstractExpr, y::AbstractExpr)
        if sign(x) == ComplexSign() || sign(y) == ComplexSign()
            error("Both the arguments should be real but these are instead $(sign(x)) and $(sign(y))")
        else
            children = (x, y)
            return new(:relative_entropy, hash(children), children, (1,1))
        end
    end
end

function sign(x::RelativeEntropyAtom)
    return NoSign()
end

function monotonicity(x::RelativeEntropyAtom)
    return (NoMonotonicity(),NoMonotonicity())
end

function curvature(x::RelativeEntropyAtom)
    return ConvexVexity()
end

function evaluate(e::RelativeEntropyAtom)
    x = vectorize(evaluate(e.children[1]))
    y = vectorize(evaluate(e.children[2]))
    if any(isnan, y) return Inf end

    out = x.*log.(x./y)
    # fix value when x=0:
    # out will only be NaN if x=0, in which case the correct value is 0
    out[isnan.(out)] .= 0
    return sum(out)
end

function conic_form!(e::RelativeEntropyAtom, unique_conic_forms::UniqueConicForms)
    if !has_conic_form(unique_conic_forms, e)
        # transform to conic form:
        # x log x/y <= z
        # x log y/x >= -z
        # log y/x >= -z/x
        # y/x >= exp(-z/x)
        # y >= x exp(-z/x)
        # and cf the standard form for the exponential cone {(x,y,z): y*exp(x/y) <= z}
        z = Variable(e.size)
        x = e.children[1]
        y = e.children[2]
        objective = conic_form!(z, unique_conic_forms)
        for i=1:size(x,1)
            for j=1:size(x,2)
                conic_form!(ExpConstraint(-z[i,j], x[i,j], y[i,j]), unique_conic_forms)
            end
        end
        # need to constrain x>=0 and y>0.
        # x>=0 we get for free from the form of the exponential cone, so just add
        conic_form!(y>=0, unique_conic_forms) # nb we don't know how to ask for strict inequality
        cache_conic_form!(unique_conic_forms, e, objective)
    end
    return get_conic_form(unique_conic_forms, e)
end


function test(optimizer)
    model = Convex.Context{Float64}(optimizer).model
    u = MOI.add_variable(model)
    v = MOI.add_variable(model)
    w = MOI.add_variable(model)
    # f = MOI.VectorOfVariables([u, v, w])
    f = MOI.VectorAffineFunction{Float64}([MOI.VectorAffineTerm{Float64}(1, MOI.ScalarAffineTerm{Float64}(1.0, u)), MOI.VectorAffineTerm{Float64}(2, MOI.ScalarAffineTerm{Float64}(1.0, v)), MOI.VectorAffineTerm{Float64}(3, MOI.ScalarAffineTerm{Float64}(1.0, w))], 
    [0.0, 0.0, 0.0])
    MOI.add_constraint(model, f, MOI.RelativeEntropyCone(1))

    g2 = MOI.VectorAffineFunction{Float64}([MOI.VectorAffineTerm{Float64}(1, MOI.ScalarAffineTerm{Float64}(1.0, v))], 
    [-5.0])
    MOI.add_constraint(model, g2, MOI.Zeros(1))

    g1 = MOI.VectorAffineFunction{Float64}([MOI.VectorAffineTerm{Float64}(1, MOI.ScalarAffineTerm{Float64}(1.0, w))], 
    [-10.0])
    MOI.add_constraint(model, g1, MOI.Nonpositives(1))
    # MOI.add_constraint(model, MOI.SingleVariable(v), MOI.EqualTo(5.0))
    # MOI.add_constraint(model, MOI.SingleVariable(w), MOI.LessThan(10.0))
    # o=MOI.SingleVariable(u)
    o = MOI.ScalarAffineFunction{Float64}(MOI.ScalarAffineTerm{Float64}[MOI.ScalarAffineTerm{Float64}(-1.0, u)], 0.0)
    MOI.set(model, MOI.ObjectiveFunction{typeof(o)}(), o)
    MOI.set(model, MOI.ObjectiveSense(), MOI.MAX_SENSE)
    MOI.optimize!(model)
    obj_val= MOI.get(model, MOI.ObjectiveValue())
    u_val = MOI.get(model, MOI.VariablePrimal(), u)
    return obj_val, u_val
end

function template(e::RelativeEntropyAtom, context::Context{T}) where {T}
    # relative_entropy(x,y) = sum_i( x_i log (x_i/y_i) ) 
    w = template(e.children[1], context)
    v = template(e.children[2], context)
    u = template(Variable(), context)
    f = operate(vcat, T, u, v, w)
    d = MOI.output_dimension(w)
    @assert d == MOI.output_dimension(v)
    MOI_add_constraint(context.model, f, MOI.RelativeEntropyCone(2d+1))
    return u
end
# fallback
operate(op::F, ::Type{T}, args...) where {F, T} = op(args...)


relative_entropy(x::AbstractExpr, y::AbstractExpr) = RelativeEntropyAtom(x, y)
# y*log(x/y)
log_perspective(x::AbstractExpr, y::AbstractExpr) = -relative_entropy(y, x)
