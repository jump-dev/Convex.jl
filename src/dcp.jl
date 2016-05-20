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

import Base.-, Base.+, Base.*
export Vexity, ConstVexity, AffineVexity, ConvexVexity, ConcaveVexity, NotDcp
export Monotonicity, Nonincreasing, Nondecreasing, NoMonotonicity
export Sign, Positive, Negative, NoSign, ComplexSign
export -, +, *

# Vexity subtypes
abstract Vexity
type ConstVexity <: Vexity              end
type AffineVexity <: Vexity             end
type ConvexVexity <: Vexity             end
type ConcaveVexity <: Vexity            end

type NotDcp <: Vexity
	function NotDcp()
		warn("Expression not DCP compliant. Trying to solve non-DCP compliant problems can lead to unexpected behavior.")
    return new()
	end
end

# Monotonocity subtypes
abstract Monotonicity
type Nonincreasing <: Monotonicity      end
type Nondecreasing <: Monotonicity      end
type ConstMonotonicity <: Monotonicity  end
type NoMonotonicity <: Monotonicity     end

# Sign subtypes
abstract Sign
type Positive <: Sign                   end
type Negative <: Sign                   end
type NoSign <: Sign                     end

# New coded

# Also create a new subtype of Sign "NotDefined to handle the ComplexSign case"
type ComplexSign <: Sign                 end

# # New code
# # Domain Subtypes
# abstract Domain
# type Real <: Domain                    end
# type ComplexSign <: Domain                 end


-(v::Vexity) = v
-(v::ConcaveVexity) = ConvexVexity()
-(v::ConvexVexity) = ConcaveVexity()

-(m::Monotonicity) = m
-(m::Nonincreasing) = Nondecreasing()
-(m::Nondecreasing) = Nonincreasing()

-(s::Sign) = s
-(s::Positive) = Negative()
-(s::Negative) = Positive()
# New code
-(s::ComplexSign) = ComplexSign()

# # New Code
# # Adding rule for Domain
# -(d::Domain) = d

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

+(s::Positive, t::Positive) = s
+(s::Negative, t::Negative) = s
+(s::Positive, t::Negative) = NoSign()
+(s::Negative, t::Positive) = NoSign()
+(s::NoSign, t::NoSign) = s
+(s::NoSign, t::Positive) = s
+(s::NoSign, t::Negative) = s
# New code
# Any sign + ComplexSign = ComplexSign
+(s::Sign, t::ComplexSign) = t
+(t::ComplexSign, s::Sign) = t

*(s::NoSign, t::NoSign) = s
*(s::NoSign, t::Positive) = s
*(s::Positive, t::NoSign) = t
*(s::NoSign, t::Negative) = s
*(s::Negative, t::NoSign) = t
*(s::Positive, t::Positive) = s
*(s::Positive, t::Negative) = t
*(s::Negative, t::Positive) = s
*(s::Negative, t::Negative) = Positive()
# New code 
# ComplexSign * Any Sign = NotDefined(Though ComplexSign and its conjugate is real but we ignore that case)
*(t::ComplexSign, s::ComplexSign) = t
*(t::ComplexSign, s::Sign) = t
*(s::Sign, t::ComplexSign) = t

*(s::Positive, m::Monotonicity) = m
*(s::Negative, m::Monotonicity) = -m
*(s::NoSign, m::Monotonicity) = NoMonotonicity()
# New Code
# ComplexSign * Any monotonivity = NoMonotonicity
*(s::ComplexSign, m::Monotonicity) = NoMonotonicity()
*(m::Monotonicity, s::Sign) = s * m

*(m::Nondecreasing, v::Vexity) = v
*(m::Nonincreasing, v::Vexity) = -v
*(m::NoMonotonicity, v::Vexity) = v
*(m::NoMonotonicity, v::ConvexVexity) = NotDcp()
*(m::NoMonotonicity, v::ConcaveVexity) = NotDcp()

# New Code
# ComplexSign * Affine = Affine
# ComplexSign * Concave = NotDcp
# ComplexSign * NotDcp = NotDcp
# ComplexSign * NotDcp = NotDcp
*(s::ComplexSign, v::ConstVexity) = v
*(s::ComplexSign, v::AffineVexity) = v
*(s::ComplexSign, v::ConvexVexity) = NotDcp()
*(s::ComplexSign, v::ConcaveVexity) = NotDcp()
#*(s::ComplexSign, v::NotDcp()) = v
*(v::Vexity, s::ComplexSign) = s*v