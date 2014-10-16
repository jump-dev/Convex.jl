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

export Vexity, ConstVexity, AffineVexity, ConvexVexity, ConcaveVexity, NotDcp
export Monotonicity, Nonincreasing, Nondecreasing, NoMonotonicity
export Sign, Positive, Negative, NoSign, Semidefinite
export -, +, *

# Vexity subtypes
abstract Vexity
type ConstVexity <: Vexity              end
type AffineVexity <: Vexity                   end
type ConvexVexity <: Vexity                   end
type ConcaveVexity <: Vexity                  end
type NotDcp <: Vexity                   end

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
type Semidefinite <: Sign               end

-(v::Vexity) = v
-(v::ConcaveVexity) = ConvexVexity()
-(v::ConvexVexity) = ConcaveVexity()

-(m::Monotonicity) = m
-(m::Nonincreasing) = Nondecreasing()
-(m::Nondecreasing) = Nonincreasing()

-(s::Sign) = s
-(s::Positive) = Negative()
-(s::Negative) = Positive()

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
+(s::NoSign, t::Sign) = s
+(s::Sign, t::NoSign) = t

*(s::NoSign, t::NoSign) = s
*(s::NoSign, t::Sign) = s
*(s::Sign, t::NoSign) = t
*(s::Positive, t::Positive) = s
*(s::Positive, t::Negative) = t
*(s::Negative, t::Positive) = s
*(s::Negative, t::Negative) = Positive()

*(s::Positive, m::Monotonicity) = m
*(s::Negative, m::Monotonicity) = -m
*(s::NoSign, m::Monotonicity) = NoMonotonicity()
*(m::Monotonicity, s::Sign) = s * m

*(m::Nondecreasing, v::Vexity) = v
*(m::Nonincreasing, v::Vexity) = -v
*(m::NoMonotonicity, v::Vexity) = v
*(m::NoMonotonicity, v::ConvexVexity) = NotDcp()
*(m::NoMonotonicity, v::ConcaveVexity) = NotDcp()
