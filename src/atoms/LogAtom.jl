# Copyright (c) 2014: Madeleine Udell and contributors
#
# Use of this source code is governed by a BSD-style license that can be found
# in the LICENSE file or at https://opensource.org/license/bsd-2-clause

mutable struct LogAtom <: AbstractExpr
    children::Tuple{AbstractExpr}
    size::Tuple{Int,Int}

    function LogAtom(x::AbstractExpr)
        if sign(x) == ComplexSign()
            error(
                "[LogAtom] the argument should be real but it's instead complex",
            )
        end
        return new((x,), x.size)
    end
end

head(io::IO, ::LogAtom) = print(io, "log")

Base.sign(::LogAtom) = NoSign()

monotonicity(::LogAtom) = (Nondecreasing(),)

curvature(::LogAtom) = ConcaveVexity()

evaluate(x::LogAtom) = log.(evaluate(x.children[1]))

Base.log(x::AbstractExpr) = LogAtom(x)

function new_conic_form!(context::Context, e::LogAtom)
    # log(z) \geq x  <=> (x,1,z) \in ExpCone
    z = e.children[1]
    m, n = size(z)
    x = Variable(m, n)
    for i in 1:m, j in 1:n
        f = vcat(x[i, j], 1, z[i, j])
        add_constraint!(context, GenericConstraint{MOI.ExponentialCone}(f))
    end
    return conic_form!(context, x)
end
