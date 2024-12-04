# Copyright (c) 2014: Madeleine Udell and contributors
#
# Use of this source code is governed by a BSD-style license that can be found
# in the LICENSE file or at https://opensource.org/license/bsd-2-clause

#############################################################################
# dcp.jl
# This file handles the basic rules on interactions of mathematical expressions
# to create new expressions.
#
# For example: negative of a concaveVexity expression is convexVexity, or multiplication
# of two positive expressions continue to be positive.
#
# See: http://dcp.stanford.edu/rules or the original paper at
# http://web.stanford.edu/~boyd/papers/disc_cvx_prog.html
#############################################################################

# This is a subtype of MOI.UnsupportedError so that non-DCP
# errors trigger the correct paths in MOI.Test for Convex.Optimizer.
struct DCPViolationError <: MOI.UnsupportedError end

function Base.showerror(io::IO, ::DCPViolationError)
    return print(
        io,
        "DCPViolationError: Expression not DCP compliant. This either means " *
        "that your problem is not convex, or that we could not prove it " *
        "was convex using the rules of disciplined convex programming. For a " *
        "list of supported operations, see https://jump.dev/Convex.jl/stable/operations/. " *
        "For help writing your problem as a disciplined convex program, " *
        "please post a reproducible example on https://jump.dev/forum.",
    )
end

abstract type Vexity end

struct ConstVexity <: Vexity end
struct AffineVexity <: Vexity end
struct ConvexVexity <: Vexity end
struct ConcaveVexity <: Vexity end
struct NotDcp <: Vexity end

Base.broadcastable(v::Vexity) = Ref(v)

abstract type Monotonicity end

struct Nonincreasing <: Monotonicity end
struct Nondecreasing <: Monotonicity end
struct ConstMonotonicity <: Monotonicity end
struct NoMonotonicity <: Monotonicity end

abstract type Sign end

struct Positive <: Sign end
struct Negative <: Sign end
struct NoSign <: Sign end
struct ComplexSign <: Sign end

Base.broadcastable(s::Sign) = Ref(s)

# -(::Vexity)

Base.:-(v::Vexity) = v
Base.:-(::ConcaveVexity) = ConvexVexity()
Base.:-(::ConvexVexity) = ConcaveVexity()

# -(::Monotonicity)

Base.:-(m::Monotonicity) = m
Base.:-(::Nonincreasing) = Nondecreasing()
Base.:-(::Nondecreasing) = Nonincreasing()

# -(::Sign)

Base.:-(s::Sign) = s
Base.:-(::Positive) = Negative()
Base.:-(::Negative) = Positive()

# +(::Sign)

Base.:+(s::Sign) = s

# +(::Vexity, ::Vexity)

Base.:+(v::ConvexVexity, ::ConvexVexity) = v
Base.:+(v::ConvexVexity, ::AffineVexity) = v
Base.:+(::ConvexVexity, ::ConcaveVexity) = NotDcp()
Base.:+(v::ConcaveVexity, ::ConcaveVexity) = v
Base.:+(v::ConcaveVexity, ::AffineVexity) = v
Base.:+(::ConcaveVexity, ::ConvexVexity) = NotDcp()
Base.:+(v::ConstVexity, ::ConstVexity) = v
Base.:+(::ConstVexity, v::NotDcp) = v
Base.:+(::ConstVexity, v::Vexity) = v
Base.:+(::AffineVexity, v::ConvexVexity) = v
Base.:+(::AffineVexity, v::ConcaveVexity) = v
Base.:+(v::AffineVexity, ::AffineVexity) = v
Base.:+(v::NotDcp, ::ConstVexity) = v
Base.:+(v::NotDcp, ::NotDcp) = v
Base.:+(v::NotDcp, ::Vexity) = v
Base.:+(v::Vexity, ::ConstVexity) = v
Base.:+(::Vexity, v::NotDcp) = v

# +(::Sign, ::Sign)

Base.:+(s::Positive, ::Positive) = s
Base.:+(s::Negative, ::Negative) = s
Base.:+(::Positive, ::Negative) = NoSign()
Base.:+(::Negative, ::Positive) = NoSign()
Base.:+(s::NoSign, ::NoSign) = s
Base.:+(s::NoSign, ::Sign) = s
Base.:+(::Sign, s::NoSign) = s
Base.:+(s::ComplexSign, ::ComplexSign) = s
Base.:+(::NoSign, s::ComplexSign) = s
Base.:+(s::ComplexSign, ::NoSign) = s
Base.:+(::Sign, s::ComplexSign) = s
Base.:+(s::ComplexSign, ::Sign) = s

# *(::Sign, ::Sign)

Base.:*(s::Positive, ::Positive) = s
Base.:*(::Negative, ::Negative) = Positive()
Base.:*(::Positive, s::Negative) = s
Base.:*(s::Negative, ::Positive) = s
Base.:*(s::NoSign, ::NoSign) = s
Base.:*(s::NoSign, ::Sign) = s
Base.:*(::Sign, s::NoSign) = s
Base.:*(s::ComplexSign, ::ComplexSign) = s
Base.:*(::NoSign, s::ComplexSign) = s
Base.:*(s::ComplexSign, ::NoSign) = s
Base.:*(::Sign, s::ComplexSign) = s
Base.:*(s::ComplexSign, ::Sign) = s

# *(::Sign, ::Monotonicity)

Base.:*(::Positive, m::Monotonicity) = m
Base.:*(::Negative, m::Monotonicity) = -m
Base.:*(::NoSign, ::Monotonicity) = NoMonotonicity()
Base.:*(::ComplexSign, ::Monotonicity) = NoMonotonicity()

# *(::Monotonicity, ::Sign)

Base.:*(m::Monotonicity, s::Sign) = s * m

# *(::Monotonicity, ::Vexity)

Base.:*(::Nondecreasing, v::Vexity) = v
Base.:*(::Nonincreasing, v::Vexity) = -v
Base.:*(::ConstMonotonicity, v::Vexity) = v
Base.:*(::NoMonotonicity, v::Vexity) = v
Base.:*(::NoMonotonicity, ::ConvexVexity) = NotDcp()
Base.:*(::NoMonotonicity, ::ConcaveVexity) = NotDcp()

# *(::Vexity, ::Monotonicity)

Base.:*(v::Vexity, m::Monotonicity) = m * v

# *(::Sign, ::Vexity)

Base.:*(::ComplexSign, v::ConcaveVexity) = NotDcp()
Base.:*(::ComplexSign, v::ConvexVexity) = NotDcp()
Base.:*(::Sign, v::Vexity) = v

# *(::Vexity, ::Sign)

Base.:*(v::Vexity, s::Sign) = s * v
