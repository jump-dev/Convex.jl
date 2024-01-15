mutable struct LogSumExpAtom <: AbstractExpr
    children::Tuple{AbstractExpr}
    size::Tuple{Int,Int}

    function LogSumExpAtom(x::AbstractExpr)
        if sign(x) == ComplexSign()
            error(
                "[LogSumExpAtom] the argument should be real but it's instead complex",
            )
        end
        return new((x,), (1, 1))
    end
end

head(io::IO, ::LogSumExpAtom) = print(io, "logsumexp")

Base.sign(::LogSumExpAtom) = NoSign()

monotonicity(::LogSumExpAtom) = (Nondecreasing(),)

curvature(::LogSumExpAtom) = ConvexVexity()

function evaluate(x::LogSumExpAtom)
    _x = evaluate(x.children[1])
    max_x = maximum(_x)
    return max_x + log(sum(exp.(_x .- max_x)))
end

logsumexp(x::AbstractExpr) = LogSumExpAtom(x)

function new_conic_form!(context::Context, e::LogSumExpAtom)
    # log(sum(exp(x))) <= t  <=>  sum(exp(x)) <= exp(t) <=> sum(exp(x - t)) <= 1
    t = Variable()
    z = sum(exp(e.children[1] - t * ones(size(e.children[1]))))
    add_constraint!(context, z <= 1)
    return conic_form!(context, t)
end

function logisticloss(e::AbstractExpr)
    if length(e) == 1
        return logsumexp([e; 0])
    end
    return sum(logsumexp([e[i]; 0]) for i in 1:length(e))
end
