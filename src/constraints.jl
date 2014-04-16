export CvxConstr, ==, >=, <=, <, >

### XXX add semidefinite constraint

type CvxConstr
  head
  lhs
  rhs
  vexity
  size
  dual_value
    canon_form
  function CvxConstr(head::Symbol,lhs::AbstractCvxExpr,rhs::AbstractCvxExpr)
    # promote sizes for zero-length dimensions, and check others match
    size = promote_size(lhs,rhs)

    # check vexity
    if head == :(==)
      if lhs.vexity in (:linear,:constant)  && rhs.vexity in (:linear,:constant)
        vexity = :linear
      else
        error("equality constraints between nonlinear expressions are not DCP compliant")
      end
    elseif head == :(<=)
      if lhs.vexity in (:linear,:constant,:convex)  && rhs.vexity in (:linear,:constant,:concave)
        vexity = :convex
      else
        error("constraint is not DCP compliant")
      end
    elseif head == :(>=)
      if rhs.vexity in (:linear,:constant,:convex)  && lhs.vexity in (:linear,:constant,:concave)
        vexity = :convex
      else
        error("constraint is not DCP compliant")
      end   
    else
      error("unrecognized comparison $head")
    end

    if lhs.head == :constant && rhs.head == :constant
      error ("TODO")
    elseif lhs.head == :constant
      if head == :(>=)
        return {:coeffs => [1], :vars => [unique_id(rhs)], :constant => lhs.value, :is_eq => false}
      else 
        return {:coeffs => [-1], :vars => [unique_id(rhs)], :constant => -(lhs.value), :is_eq => (head == :(==))}
      end
    elseif rhs.head == :constant
      if head == :(>=)
        return {:coeffs => [-1], :vars => [unique_id(lhs)], :constant => -(rhs.value), :is_eq => false}
      else
        return {:coeffs => [1], :vars => [unique_id(lhs)], :constant => rhs.value, :is_eq => (head == :(==))}
      end
    else
      if head == :(>=)
        return {:coeffs => [1 -1], :vars => [unique_id(rhs); unique_id(lhs)], :constant => 0, :is_eq => false}
      else
        return {:coeffs => [1 -1], :vars => [unique_id(lhs); unique_id(rhs)], :constant => 0, :is_eq => (head == :(==))}
      end
    end

    return new(head,lhs,rhs,vexity,size,nothing)
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