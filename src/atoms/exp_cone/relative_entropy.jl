mutable struct RelativeEntropyAtom <: AbstractExpr
    children::Tuple{AbstractExpr,AbstractExpr}
    size::Tuple{Int,Int}

    function RelativeEntropyAtom(x::AbstractExpr, y::AbstractExpr)
        if sign(x) == ComplexSign() || sign(y) == ComplexSign()
            error(
                "Both the arguments should be real but these are instead $(sign(x)) and $(sign(y))",
            )
        end
        return new((x, y), (1, 1))
    end
end

head(io::IO, ::RelativeEntropyAtom) = print(io, "relative_entropy")

Base.sign(::RelativeEntropyAtom) = NoSign()

monotonicity(::RelativeEntropyAtom) = (NoMonotonicity(), NoMonotonicity())

curvature(::RelativeEntropyAtom) = ConvexVexity()

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

function new_conic_form!(context::Context{T}, e::RelativeEntropyAtom) where {T}
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

relative_entropy(x::AbstractExpr, y::AbstractExpr) = RelativeEntropyAtom(x, y)

# y*log(x/y)
log_perspective(x::AbstractExpr, y::AbstractExpr) = -relative_entropy(y, x)
