mutable struct AbsAtom <: AbstractExpr
    children::Tuple{AbstractExpr}
    size::Tuple{Int,Int}

    AbsAtom(x::AbstractExpr) = new((x,), x.size)
end

head(io::IO, ::AbsAtom) = print(io, "abs")

Base.sign(::AbsAtom) = Positive()

monotonicity(x::AbsAtom) = (Nondecreasing() * sign(x.children[1]),)

curvature(::AbsAtom) = ConvexVexity()

evaluate(x::AbsAtom) = abs.(evaluate(x.children[1]))

function new_conic_form!(context::Context, A::AbsAtom)
    x = only(A.children)
    t = Variable(size(x))
    t_obj = conic_form!(context, t)
    if sign(x) == ComplexSign()
        for i in 1:length(vec(t))
            add_constraint!(
                context,
                t[i] >= LinearAlgebra.norm2([real(x[i]); imag(x[i])]),
            )
        end
    else
        add_constraint!(context, t >= x)
        add_constraint!(context, t >= -x)
    end
    return t_obj
end

Base.abs(x::AbstractExpr) = AbsAtom(x)

Base.abs2(x::AbstractExpr) = square(abs(x))
