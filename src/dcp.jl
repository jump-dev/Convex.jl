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

import Base.-, Base.+, Base.*, Base./
export Vexity, ConstVexity, AffineVexity, ConvexVexity, ConcaveVexity, NotDcp
export Monotonicity, Nonincreasing, Nondecreasing, NoMonotonicity
export Sign, Positive, Negative, NoSign, ComplexSign
export -, +, *

# Vexity subtypes
abstract type Vexity end
struct ConstVexity <: Vexity              end
struct AffineVexity <: Vexity             end
struct ConvexVexity <: Vexity             end
struct ConcaveVexity <: Vexity            end

struct NotDcp <: Vexity
    function NotDcp()
        @warn """
              Expression not DCP compliant.
              Trying to solve non-DCP compliant problems can lead to unexpected behavior.
              """
        return new()
    end
end

# Monotonocity subtypes
abstract type Monotonicity end
struct Nonincreasing <: Monotonicity      end
struct Nondecreasing <: Monotonicity      end
struct ConstMonotonicity <: Monotonicity  end
struct NoMonotonicity <: Monotonicity     end

# Sign subtypes
abstract type Sign end
struct Positive <: Sign                   end
struct Negative <: Sign                   end
struct NoSign <: Sign                     end

# New coded

# Also create a new subtype of Sign "NotDefined to handle the ComplexSign case"
struct ComplexSign <: Sign                end


-(v::Vexity) = v
-(v::ConcaveVexity) = ConvexVexity()
-(v::ConvexVexity) = ConcaveVexity()

-(m::Monotonicity) = m
-(m::Nonincreasing) = Nondecreasing()
-(m::Nondecreasing) = Nonincreasing()

-(s::Sign) = s
-(s::Positive) = Negative()
-(s::Negative) = Positive()
-(s::ComplexSign) = ComplexSign()



+(v::NotDcp, w::NotDcp) = v
+(v::NotDcp, w::Vexity) = v
+(v::Vexity, w::NotDcp) = w

+(v::ConstVexity, w::ConstVexity) = v
+(v::ConstVexity, w::NotDcp) = w
+(v::NotDcp, w::ConstVexity) = v
+(v::ConstVexity, w::Vexity) = w
+(v::Vexity, w::ConstVexity) = v

+(v::AffineVexity, w::AffineVexity) = v
+(v::AffineVexity, w::ConvexVexity) = w
+(v::ConvexVexity, w::AffineVexity) = v
+(v::AffineVexity, w::ConcaveVexity) = w
+(v::ConcaveVexity, w::AffineVexity) = v

+(v::ConvexVexity, w::ConvexVexity) = v
+(v::ConcaveVexity, w::ConcaveVexity) = v
+(v::ConcaveVexity, w::ConvexVexity) = NotDcp()
+(v::ConvexVexity, w::ConcaveVexity) = NotDcp()

#+(::Convex.Positive, ::Convex.NoSign)
+(s::Positive, t::Positive) = s
+(s::Negative, t::Negative) = s
+(s::Positive, t::Negative) = NoSign()
+(s::Negative, t::Positive) = NoSign()
+(s::NoSign, t::NoSign) = s
+(s::NoSign, t::Positive) = s
+(t::Positive, s::NoSign) = s+t
+(s::NoSign, t::Negative) = s
+(t::Negative, s::NoSign) = s+t

# Any sign + ComplexSign = ComplexSign
+(s::ComplexSign, t::ComplexSign) = s
+(s::Sign, t::ComplexSign) = t
+(t::ComplexSign, s::Sign) = s+t

*(s::NoSign, t::NoSign) = s
*(s::NoSign, t::Positive) = s
*(s::Positive, t::NoSign) = t
*(s::NoSign, t::Negative) = s
*(s::Negative, t::NoSign) = t
*(s::Positive, t::Positive) = s
*(s::Positive, t::Negative) = t
*(s::Negative, t::Positive) = s
*(s::Negative, t::Negative) = Positive()

# ComplexSign * Any Sign = NotDefined(Though ComplexSign and its conjugate is real but we ignore that case)
*(t::ComplexSign, s::ComplexSign) = t
*(t::ComplexSign, s::Sign) = t
*(s::Sign, t::ComplexSign) = t

*(s::Positive, m::Monotonicity) = m
*(s::Negative, m::Monotonicity) = -m
*(s::NoSign, m::Monotonicity) = NoMonotonicity()

# ComplexSign * Any monotonivity = NoMonotonicity
*(s::ComplexSign, m::Monotonicity) = NoMonotonicity()
*(m::Monotonicity, s::Sign) = s * m

*(m::Nondecreasing, v::Vexity) = v
*(m::Nonincreasing, v::Vexity) = -v
*(m::NoMonotonicity, v::Vexity) = v
*(m::NoMonotonicity, v::ConvexVexity) = NotDcp()
*(m::NoMonotonicity, v::ConcaveVexity) = NotDcp()


# ComplexSign * Affine = Affine
# ComplexSign * Concave = NotDcp
# ComplexSign * NotDcp = NotDcp
# ComplexSign * NotDcp = NotDcp
*(s::ComplexSign, v::ConstVexity) = v
*(s::ComplexSign, v::AffineVexity) = v
*(s::ComplexSign, v::ConvexVexity) = NotDcp()
*(s::ComplexSign, v::ConcaveVexity) = NotDcp()
*(s::ComplexSign, v::NotDcp) = v
*(v::Vexity, s::ComplexSign) = s*v
