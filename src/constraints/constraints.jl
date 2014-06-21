export CanonicalConstr, CvxConstr, ==, >=, <=, <, >, +
export .==, .>=, .<=, .>, .<

# TODO: Allow matrix-vector comparisons if one dimension matches

# Holds a constraint in canonical_form
# We have an array of `coeffs` where these coefficients corresponds to `vars`
# is_eq: the constraint is either == or <=
# is_conic: specifies whether or not is conic
type CanonicalConstr
  coeffs::VecOrMatOrSparse
  vars::Array{Int64, 1}
  constant::Value
  is_eq::Bool
  is_conic::Bool
  uid::Int64

  function CanonicalConstr(coeffs::VecOrMat, vars::Array{Int64, 1}, constant::Value, is_eq::Bool, is_conic::Bool)
    this = new(coeffs, vars, constant, is_eq, is_conic)
    this.uid = unique_id(this)
    return this
  end

  function CanonicalConstr(coeffs::Number, vars::Array{Int64, 1}, constant::Value, is_eq::Bool, is_conic::Bool)
    this = new([coeffs], vars, constant, is_eq, is_conic)
    this.uid = unique_id(this)
    return this
  end

  function CanonicalConstr(coeffs::VecOrMat, vars::Int64, constant::Value, is_eq::Bool, is_conic::Bool)
    this = new(coeffs, [vars], constant, is_eq, is_conic)
    this.uid = unique_id(this)
    return this
  end

  function CanonicalConstr(coeffs::Number, vars::Int64, constant::Value, is_eq::Bool, is_conic::Bool)
    this = new([coeffs], [vars], constant, is_eq, is_conic)
    this.uid = unique_id(this)
    return this
  end
end

# A constraint between two `AbstractCvxExpr`
type CvxConstr
  head::Symbol
  lhs::AbstractCvxExpr
  rhs::AbstractCvxExpr
  vexity::Symbol
  size::(Int64, Int64)
  canon_form::Function
  canon_uid::Int64
  dual_value::ValueOrNothing

  function CvxConstr(head::Symbol, lhs::AbstractCvxExpr, rhs::AbstractCvxExpr)
    # Check vexity
    if head == :(==)
      if lhs.vexity in (:affine, :constant)  && rhs.vexity in (:affine, :constant)
        vexity = :affine
      else
        error("Equality constraints between nonaffine expressions are not DCP compliant")
      end
    elseif head == :(<=)
      if lhs.vexity in (:affine, :constant, :convex) && rhs.vexity in (:affine, :constant, :concave)
        vexity = :convex
      else
        error("Constraint is not DCP compliant")
      end
    elseif head == :(>=)
      error(">= should have been transformed to <=")
    else
      error("Unrecognized comparison $head")
    end

    # Size checks will be done in canon_form
    if lhs.size == (1, 1)
      sz = rhs.size
    else
      sz = lhs.size
    end
    constraint = new(head, lhs, rhs, vexity, sz)

    canon_form = ()->
      begin
        # If both lhs and rhs are constants, we check if the constraint is feasible
        # If it is not feasible, we throw an error. If it is feasible, the
        # canonical_form will be an empty array
        if lhs.vexity == :constant && rhs.vexity == :constant
          if head == :(==)
            if lhs.value != rhs.value
              error("Infeasible problem. $(lhs.value) != $(rhs.value)")
            end
          elseif head == :(<=)
            if !(lhs.value <= rhs.value)
              error("Infeasible problem. $(lhs.value) is not <= $(rhs.value)")
            end
          end
          return CanonicalConstr[]
        # In this case, we have lhs <= constant
        # We promote the size of the rhs if it is 1 x 1
        # If the lhs is 1 x 1, we write this as ones(...) * lhs <= constant
        # to make comparisons of the same size
        elseif rhs.vexity == :constant
          if rhs.size == (1, 1) && lhs.size != (1, 1)
            rhs = Constant(rhs.value * ones(lhs.size...), rhs.sign)
            coeffs = VecOrMatOrSparse[speye(get_vectorized_size(lhs))]
          elseif rhs.size != (1, 1) && lhs.size == (1, 1)
            coeffs = VecOrMatOrSparse[ones(get_vectorized_size(rhs), 1)]
          elseif lhs.size != rhs.size
            error("Can't compare expressions of size $(x.size) and $(y.size)")
          else
            coeffs = VecOrMatOrSparse[speye(get_vectorized_size(lhs))]
          end

          constant = vec(rhs.value)
          canon_constr = CanonicalConstr(coeffs, unique_id(lhs), constant, (head == :(==)), false)
          constraint.canon_uid = canon_constr.uid

          canon_constr_array = lhs.canon_form()
          push!(canon_constr_array, canon_constr)

        # If we have lhs <= rhs, the canonical form is [I -I] * [lhs rhs]' <= 0
        # Again, we multiply by ones(...) if the size of lhs or rhs is (1, 1)
        else
          if lhs.size == (1, 1) && rhs.size != (1, 1)
            sz = get_vectorized_size(rhs.size)
            coeffs = VecOrMatOrSparse[ones(sz, 1), -speye(sz)]
          elseif lhs.size != (1, 1) && rhs.size == (1, 1)
            sz = get_vectorized_size(lhs.size)
            coeffs = VecOrMatOrSparse[speye(sz), -ones(sz, 1)]
          elseif lhs.size != rhs.size
            error("Can't compare expressions of size $(x.size) and $(y.size)")
          else
            sz = get_vectorized_size(rhs.size)
            coeffs = VecOrMatOrSparse[speye(sz), -speye(sz)]
          end
          vars = [unique_id(lhs); unique_id(rhs)]
          constant = zeros(sz)

          canon_constr = CanonicalConstr(coeffs, vars, constant, (head == :(==)), false)
          constraint.canon_uid = canon_constr.uid

          canon_constr_array = lhs.canon_form()
          append!(canon_constr_array, rhs.canon_form())
          push!(canon_constr_array, canon_constr)
        end
        return canon_constr_array
      end

    constraint.canon_form = canon_form
    return constraint
  end
end

==(x::AbstractCvxExpr, y::AbstractCvxExpr) = CvxConstr(:(==), x, y)
>=(x::AbstractCvxExpr, y::AbstractCvxExpr) = CvxConstr(:(<=), y, x)
<=(x::AbstractCvxExpr, y::AbstractCvxExpr) = CvxConstr(:(<=), x, y)
>(x::AbstractCvxExpr, y::AbstractCvxExpr) = >=(x, y)
<(x::AbstractCvxExpr, y::AbstractCvxExpr) = <=(x, y)

==(x::Constant, y::AbstractCvxExpr) = CvxConstr(:(==), y, x)
<=(x::Constant, y::AbstractCvxExpr) = CvxConstr(:(<=), -y, -x)
>=(x::Constant, y::AbstractCvxExpr) = CvxConstr(:(<=), y, x)
>(x::Constant, y::AbstractCvxExpr) = <=(y, x)
<(x::Constant, y::AbstractCvxExpr) = >=(y, x)

==(x::Value, y::AbstractCvxExpr) = CvxConstr(:(==), y, convert(CvxExpr, x))
>=(x::Value, y::AbstractCvxExpr) = CvxConstr(:(<=), y, convert(CvxExpr, x))
<=(x::Value, y::AbstractCvxExpr) = CvxConstr(:(<=), -y, -convert(CvxExpr, x))
>(x::Value, y::AbstractCvxExpr) = <=(y, x)
<(x::Value, y::AbstractCvxExpr) = >=(y, x)

==(x::AbstractCvxExpr, y::Value)= CvxConstr(:(==), x, convert(CvxExpr, y))
>=(x::AbstractCvxExpr, y::Value) = CvxConstr(:(<=), -x, -convert(CvxExpr, y))
<=(x::AbstractCvxExpr, y::Value) = CvxConstr(:(<=), x, convert(CvxExpr, y))
>(x::AbstractCvxExpr, y::Value) = >=(x, y)
<(x::AbstractCvxExpr, y::Value) = <=(x, y)

+(constraints::Array{CvxConstr}, new_constraint::CvxConstr) =
  push!(constraints, new_constraint)
+(constraints::Array{CvxConstr}, new_constraints::Array{CvxConstr}) =
  append!(constraints, new_constraints)

# Using .<= or <= have the same behavior as of now
.==(x::AbstractCvxExpr, y::AbstractCvxExpr) = CvxConstr(:(==), x, y)
.>=(x::AbstractCvxExpr, y::AbstractCvxExpr) = CvxConstr(:(<=), y, x)
.<=(x::AbstractCvxExpr, y::AbstractCvxExpr) = CvxConstr(:(<=), x, y)
.>(x::AbstractCvxExpr, y::AbstractCvxExpr) = >=(x, y)
.<(x::AbstractCvxExpr, y::AbstractCvxExpr) = <=(x, y)

.==(x::Constant, y::AbstractCvxExpr) = CvxConstr(:(==), y, x)
.<=(x::Constant, y::AbstractCvxExpr) = CvxConstr(:(<=), -y, -x)
.>=(x::Constant, y::AbstractCvxExpr) = CvxConstr(:(<=), y, x)
.>(x::Constant, y::AbstractCvxExpr) = <=(y, x)
.<(x::Constant, y::AbstractCvxExpr) = >=(y, x)

# The following needs AbstractArray to be handled separately due to ambiguity
# that otherwise arises with Julia
.==(x::AbstractArray, y::AbstractCvxExpr) = CvxConstr(:(==), y, convert(CvxExpr, x))
.==(x::Value, y::AbstractCvxExpr) = CvxConstr(:(==), y, convert(CvxExpr, x))

.>=(x::Value, y::AbstractCvxExpr) = CvxConstr(:(<=), y, convert(CvxExpr, x))

.<=(x::AbstractArray, y::AbstractCvxExpr) = CvxConstr(:(<=), -y, -convert(CvxExpr, x))
.<=(x::Value, y::AbstractCvxExpr) = CvxConstr(:(<=), -y, -convert(CvxExpr, x))

.>(x::Value, y::AbstractCvxExpr) = <=(y, x)

.<(x::AbstractArray, y::AbstractCvxExpr) = >=(y, x)
.<(x::Value, y::AbstractCvxExpr) = >=(y, x)

.==(x::AbstractCvxExpr, y::AbstractArray)= CvxConstr(:(==), x, convert(CvxExpr, y))
.==(x::AbstractCvxExpr, y::Value)= CvxConstr(:(==), x, convert(CvxExpr, y))

.>=(x::AbstractCvxExpr, y::Value) = CvxConstr(:(<=), -x, -convert(CvxExpr, y))

.<=(x::AbstractCvxExpr, y::AbstractArray) = CvxConstr(:(<=), x, convert(CvxExpr, y))
.<=(x::AbstractCvxExpr, y::Value) = CvxConstr(:(<=), x, convert(CvxExpr, y))

.>(x::AbstractCvxExpr, y::Value) = >=(x, y)

.<(x::AbstractCvxExpr, y::AbstractArray) = <=(x, y)
.<(x::AbstractCvxExpr, y::Value) = <=(x, y)
