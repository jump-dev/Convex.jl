mutable struct EntropyAtom <: AbstractExpr
    children::Tuple{AbstractExpr}
    size::Tuple{Int,Int}

    function EntropyAtom(x::AbstractExpr)
        if sign(x) == ComplexSign()
            error(
                "[EntropyAtom] the argument should be real but it's instead complex",
            )
        end
        # TODO check positivity or enforce it
        return new((x,), size(x))
    end
end

head(io::IO, ::EntropyAtom) = print(io, "entropy")

Base.sign(::EntropyAtom) = NoSign()

monotonicity(::EntropyAtom) = (NoMonotonicity(),)

curvature(::EntropyAtom) = ConcaveVexity()

function evaluate(x::EntropyAtom)
    c = evaluate(x.children[1])
    return -c .* log.(c)
end

entropy(x::AbstractExpr) = sum(EntropyAtom(x))

entropy_elementwise(x::AbstractExpr) = EntropyAtom(x)

function new_conic_form!(context::Context, e::EntropyAtom)
    # -x log x >= t  <=>  x exp(t/x) <= 1  <==>  (t,x,1) in exp cone
    t = Variable(e.size)
    x = e.children[1]
    add_constraint!(
        context,
        GenericConstraint{MOI.ExponentialCone}(vcat(t, x, 1)),
    )
    return conic_form!(context, t)
end
