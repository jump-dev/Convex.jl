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

# Vexity subtypes
abstract type Vexity end
Base.broadcastable(v::Vexity) = Ref(v)
struct ConstVexity <: Vexity end
struct AffineVexity <: Vexity end
struct ConvexVexity <: Vexity end
struct ConcaveVexity <: Vexity end

struct DCPViolationError <: Exception end

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

struct NotDcp <: Vexity end

# Monotonocity subtypes
abstract type Monotonicity end
struct Nonincreasing <: Monotonicity end
struct Nondecreasing <: Monotonicity end
struct ConstMonotonicity <: Monotonicity end
struct NoMonotonicity <: Monotonicity end

# Sign subtypes
abstract type Sign end
Base.broadcastable(s::Sign) = Ref(s)
struct Positive <: Sign end
struct Negative <: Sign end
struct NoSign <: Sign end

# New coded

# Also create a new subtype of Sign "NotDefined to handle the ComplexSign case"
struct ComplexSign <: Sign end

Base.:-(v::Vexity) = v
Base.:-(v::ConcaveVexity) = ConvexVexity()
Base.:-(v::ConvexVexity) = ConcaveVexity()

Base.:-(m::Monotonicity) = m
Base.:-(m::Nonincreasing) = Nondecreasing()
Base.:-(m::Nondecreasing) = Nonincreasing()

Base.:-(s::Sign) = s
Base.:-(s::Positive) = Negative()
Base.:-(s::Negative) = Positive()
Base.:-(s::ComplexSign) = ComplexSign()

Base.:+(v::NotDcp, w::NotDcp) = v
Base.:+(v::NotDcp, w::Vexity) = v
Base.:+(v::Vexity, w::NotDcp) = w

Base.:+(v::ConstVexity, w::ConstVexity) = v
Base.:+(v::ConstVexity, w::NotDcp) = w
Base.:+(v::NotDcp, w::ConstVexity) = v
Base.:+(v::ConstVexity, w::Vexity) = w
Base.:+(v::Vexity, w::ConstVexity) = v

Base.:+(v::AffineVexity, w::AffineVexity) = v
Base.:+(v::AffineVexity, w::ConvexVexity) = w
Base.:+(v::ConvexVexity, w::AffineVexity) = v
Base.:+(v::AffineVexity, w::ConcaveVexity) = w
Base.:+(v::ConcaveVexity, w::AffineVexity) = v

Base.:+(v::ConvexVexity, w::ConvexVexity) = v
Base.:+(v::ConcaveVexity, w::ConcaveVexity) = v
Base.:+(v::ConcaveVexity, w::ConvexVexity) = NotDcp()
Base.:+(v::ConvexVexity, w::ConcaveVexity) = NotDcp()

#Base.:+(::Convex.Positive, ::Convex.NoSign)
Base.:+(s::Sign) = s
Base.:+(s::Positive, t::Positive) = s
Base.:+(s::Negative, t::Negative) = s
Base.:+(s::Positive, t::Negative) = NoSign()
Base.:+(s::Negative, t::Positive) = NoSign()
Base.:+(s::NoSign, t::NoSign) = s
Base.:+(s::NoSign, t::Positive) = s
Base.:+(t::Positive, s::NoSign) = s + t
Base.:+(s::NoSign, t::Negative) = s
Base.:+(t::Negative, s::NoSign) = s + t

# Any sign + ComplexSign = ComplexSign
Base.:+(s::ComplexSign) = s
Base.:+(s::ComplexSign, t::ComplexSign) = s
Base.:+(s::Sign, t::ComplexSign) = t
Base.:+(t::ComplexSign, s::Sign) = s + t

Base.:*(s::NoSign, t::NoSign) = s
Base.:*(s::NoSign, t::Positive) = s
Base.:*(s::Positive, t::NoSign) = t
Base.:*(s::NoSign, t::Negative) = s
Base.:*(s::Negative, t::NoSign) = t
Base.:*(s::Positive, t::Positive) = s
Base.:*(s::Positive, t::Negative) = t
Base.:*(s::Negative, t::Positive) = s
Base.:*(s::Negative, t::Negative) = Positive()

# ComplexSign * Any Sign = NotDefined(Though ComplexSign and its conjugate is real but we ignore that case)
Base.:*(t::ComplexSign, s::ComplexSign) = t
Base.:*(t::ComplexSign, s::Sign) = t
Base.:*(s::Sign, t::ComplexSign) = t

Base.:*(s::Positive, m::Monotonicity) = m
Base.:*(s::Negative, m::Monotonicity) = -m
Base.:*(s::NoSign, m::Monotonicity) = NoMonotonicity()

# ComplexSign * Any monotonivity = NoMonotonicity
Base.:*(s::ComplexSign, m::Monotonicity) = NoMonotonicity()
Base.:*(m::Monotonicity, s::Sign) = s * m

Base.:*(m::Nondecreasing, v::Vexity) = v
Base.:*(m::Nonincreasing, v::Vexity) = -v
Base.:*(m::NoMonotonicity, v::Vexity) = v
Base.:*(m::NoMonotonicity, v::ConvexVexity) = NotDcp()
Base.:*(m::NoMonotonicity, v::ConcaveVexity) = NotDcp()

# ComplexSign * Affine = Affine
# ComplexSign * Concave = NotDcp
# ComplexSign * NotDcp = NotDcp
# ComplexSign * NotDcp = NotDcp
Base.:*(s::ComplexSign, v::ConstVexity) = v
Base.:*(s::ComplexSign, v::AffineVexity) = v
Base.:*(s::ComplexSign, v::ConvexVexity) = NotDcp()
Base.:*(s::ComplexSign, v::ConcaveVexity) = NotDcp()
Base.:*(s::ComplexSign, v::NotDcp) = v
Base.:*(v::Vexity, s::ComplexSign) = s * v
