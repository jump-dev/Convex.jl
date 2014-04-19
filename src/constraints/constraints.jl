export CvxConstr, ==, >=, <=, <, >

# TODO: Break down constraints.jl into multiple files- for equality/exponential,
# SOCP, SDP etc constraints

# TODO: CVX constraint should be an abstract class and children should be stuff
# like CVXEqualityConstr. Read:
# http://docs.julialang.org/en/release-0.2/manual/performance-tips/#break-functions-into-multiple-definitions
type CvxConstr
  head
  lhs
  rhs
  vexity
  size
  dual_value
  canon_form::Function
  function CvxConstr(head::Symbol, lhs::AbstractCvxExpr, rhs::AbstractCvxExpr)
    # promote sizes for zero-length dimensions, and check others match
    size = promote_size(lhs, rhs)

    # check vexity
    if head == :(==)
      if lhs.vexity in (:linear, :constant)  && rhs.vexity in (:linear, :constant)
        vexity = :linear
      else
        error("equality constraints between nonlinear expressions are not DCP compliant")
      end
    elseif head == :(<=)
      if lhs.vexity in (:linear, :constant, :convex) && rhs.vexity in (:linear, :constant, :concave)
        vexity = :convex
      else
        error("constraint is not DCP compliant")
      end
    elseif head == :(>=)
      if rhs.vexity in (:linear, :constant, :convex)  && lhs.vexity in (:linear, :constant, :concave)
        vexity = :convex
      else
        error("constraint is not DCP compliant")
      end
    else
      error("unrecognized comparison $head")
    end

    canon_form = ()->
      begin
        if lhs.head == :constant && rhs.head == :constant
          error ("TODO")
        elseif lhs.head == :constant
          if head == :(>=)
            return {
              :coeffs => [1],
              :vars => [unique_id(rhs)],
              :constant => lhs.value,
              :is_eq => false,
              :size => rhs.size
            }
          else
            return {
              :coeffs => [-1],
              :vars => [unique_id(rhs)],
              :constant => -(lhs.value),
              :is_eq => (head == :(==)),
              :size => rhs.size
            }
          end
        elseif rhs.head == :constant
          if head == :(>=)
            return {
              :coeffs => [-1],
              :vars => [unique_id(lhs)],
              :constant => -(rhs.value),
              :is_eq => false,
              :size => lhs.size
            }
          else
            return {
              :coeffs => [1],
              :vars => [unique_id(lhs)],
              :constant => rhs.value,
              :is_eq => (head == :(==)),
              :size => lhs.size
            }
          end
        else
          if head == :(>=)
            # TODO: fix size
            return {
              :coeffs => [1 -1],
              :vars => [unique_id(rhs); unique_id(lhs)],
              :constant => 0,
              :is_eq => false,
              :size => maximum(lhs.size)
            }
          else
            return {
              :coeffs => [1 -1],
              :vars => [unique_id(lhs); unique_id(rhs)],
              :constant => 0,
              :is_eq => (head == :(==)),
              :size => maximum(lhs.size)
            }
          end
        end
      end

    return new(head,lhs,rhs,vexity,size,nothing,canon_form)
  end
end

==(x::AbstractCvxExpr,y::AbstractCvxExpr) = CvxConstr(:(==),x,y)
>=(x::AbstractCvxExpr,y::AbstractCvxExpr) = CvxConstr(:(>=),x,y)
<=(x::AbstractCvxExpr,y::AbstractCvxExpr) = CvxConstr(:(<=),x,y)
>(x::AbstractCvxExpr,y::AbstractCvxExpr) = CvxConstr(:(>=),x,y)
<(x::AbstractCvxExpr,y::AbstractCvxExpr) = CvxConstr(:(<=),x,y)
==(x::Value,y::AbstractCvxExpr) = CvxConstr(:(==),y,convert(CvxExpr,x))
>=(x::Value,y::AbstractCvxExpr) = CvxConstr(:(<=),y,convert(CvxExpr,x))
<=(x::Value,y::AbstractCvxExpr) = CvxConstr(:(>=),y,convert(CvxExpr,x))
>(x::Value,y::AbstractCvxExpr) = CvxConstr(:(<=),y,convert(CvxExpr,x))
<(x::Value,y::AbstractCvxExpr) = CvxConstr(:(>=),y,convert(CvxExpr,x))
==(x::AbstractCvxExpr,y::Value)= CvxConstr(:(==),x,convert(CvxExpr,y))
>=(x::AbstractCvxExpr,y::Value) = CvxConstr(:(>=),x,convert(CvxExpr,y))
<=(x::AbstractCvxExpr,y::Value) = CvxConstr(:(<=),x,convert(CvxExpr,y))
>(x::AbstractCvxExpr,y::Value) = CvxConstr(:(>=),x,convert(CvxExpr,y))
<(x::AbstractCvxExpr,y::Value) = CvxConstr(:(<=),x,convert(CvxExpr,y))
