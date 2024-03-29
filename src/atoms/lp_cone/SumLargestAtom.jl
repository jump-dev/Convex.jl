mutable struct SumLargestAtom <: AbstractExpr
    children::Tuple{AbstractExpr}
    size::Tuple{Int,Int}
    k::Int

    function SumLargestAtom(x::AbstractExpr, k::Int)
        if sign(x) == ComplexSign()
            error(
                "[SumLargestAtom] argument should be real instead it is $(sign(x))",
            )
        elseif k <= 0
            error(
                "[SumLargestAtom] sumlargest and sumsmallest only support positive values of k",
            )
        elseif k > length(x)
            error(
                "[SumLargestAtom] k cannot be larger than the number of entries in x",
            )
        end
        return new((x,), (1, 1), k)
    end
end

head(io::IO, ::SumLargestAtom) = print(io, "sumlargest")

Base.sign(x::SumLargestAtom) = sign(x.children[1])

monotonicity(::SumLargestAtom) = (Nondecreasing(),)

curvature(::SumLargestAtom) = ConvexVexity()

function evaluate(x::SumLargestAtom)
    return sum(sort(vec(evaluate(x.children[1])); rev = true)[1:x.k])
end

function new_conic_form!(context::Context, x::SumLargestAtom)
    c = x.children[1]
    t = Variable(size(c))
    q = Variable()
    # sum k largest given by the solution to
    # minimize sum(t) + k*q
    # subject to c <= t + q, t >= 0
    add_constraint!(context, c <= t + q)
    add_constraint!(context, t >= 0)
    return conic_form!(context, sum(t) + x.k * q)
end

function sumlargest(x::AbstractExpr, k::Int)
    return k == 0 ? Constant(0) : SumLargestAtom(x, k)
end

function sumsmallest(x::AbstractExpr, k::Int)
    return k == 0 ? Constant(0) : -SumLargestAtom(-x, k)
end
