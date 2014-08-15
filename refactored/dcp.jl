export Vexity, ConstVexity, Affine, Convex, Concave, NoVexity
export Monotonicity, Nonincreasing, Nondecreasing, NoMonotonicity
export Sign, Positive, Negative, NoSign
export -, +, *

abstract Vexity
abstract Monotonicity
abstract Sign


type ConstVexity <: Vexity
end
type Affine <: Vexity
end
type Convex <: Vexity
end
type Concave <: Vexity
end
type NoVexity <: Vexity
end

type Nonincreasing <: Monotonicity
end
type Nondecreasing <: Monotonicity
end
type ConstMonotonicity <: Monotonicity
end
type NoMonotonicity <: Monotonicity
end

type Positive <: Sign
end
type Negative <: Sign
end
type NoSign <: Sign
end

-(v::Vexity) = v
-(v::Concave) = Convex()
-(v::Convex) = Concave()

-(m::Monotonicity) = m
-(m::Nonincreasing) = Nondecreasing()
-(m::Nondecreasing) = Nonincreasing()

-(s::Sign) = s
-(s::Positive) = Negative()
-(s::Negative) = Positive()

+(v::NoVexity, w::NoVexity) = v
+(v::NoVexity, w::Vexity) = v
+(v::Vexity, w::NoVexity) = w

+(v::ConstVexity, w::ConstVexity) = v
+(v::ConstVexity, w::NoVexity) = w
+(v::NoVexity, w::ConstVexity) = v
+(v::ConstVexity, w::Vexity) = w
+(v::Vexity, w::ConstVexity) = v

+(v::Affine, w::Affine) = v
+(v::Affine, w::Convex) = w
+(v::Convex, w::Affine) = v
+(v::Affine, w::Concave) = w
+(v::Concave, w::Affine) = v

+(v::Convex, w::Convex) = v
+(v::Concave, w::Concave) = v
+(v::Concave, w::Convex) = NoVexity()
+(v::Convex, w::Concave) = NoVexity()


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

*(m::Nondecreasing, v::Vexity) = v
*(m::Nonincreasing, v::Vexity) = -v
*(m::NoMonotonicity, v::Vexity) = v
*(m::NoMonotonicity, v::Convex) = NoVexity()
*(m::NoMonotonicity, v::Concave) = NoVexity()
