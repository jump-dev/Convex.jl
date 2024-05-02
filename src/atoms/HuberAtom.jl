# Copyright (c) 2014: Madeleine Udell and contributors
#
# Use of this source code is governed by a BSD-style license that can be found
# in the LICENSE file or at https://opensource.org/license/bsd-2-clause

mutable struct HuberAtom <: AbstractExpr
    children::Tuple{AbstractExpr}
    size::Tuple{Int,Int}
    M::Real

    function HuberAtom(x::AbstractExpr, M::Real)
        if sign(x) == ComplexSign()
            error("[HuberAtom] argument must be real")
        elseif M <= 0
            error("[HuberAtom] parameter must by a positive scalar. Got `M=$M`")
        end
        return new((x,), x.size, M)
    end
end

head(io::IO, ::HuberAtom) = print(io, "huber")

Base.sign(::HuberAtom) = Positive()

monotonicity(x::HuberAtom) = (Nondecreasing() * sign(x.children[1]),)

curvature(::HuberAtom) = ConvexVexity()

function evaluate(x::HuberAtom)
    c = evaluate(x.children[1])
    if c isa Number
        c = [c]
    end
    for i in 1:length(c)
        if c[i] <= x.M
            c[i] = c[i]^2
        else
            c[i] = 2 * x.M * c[i] - x.M^2
        end
    end
    return c
end

function new_conic_form!(context::Context, x::HuberAtom)
    c = x.children[1]
    s = Variable(c.size)
    n = Variable(c.size)
    add_constraint!(context, c == s + n)
    return conic_form!(context, square(s) + 2 * x.M * abs(n))
end
