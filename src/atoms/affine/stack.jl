import Base.vcat, Base.hcat
export vcat, hcat, VcatAtom, HcatAtom
export sign, curvature, monotonicity, evaluate, conic_form!

type HcatAtom <: AbstractExpr
  head::Symbol
  id_hash::UInt64
  children::Tuple
  size::Tuple{Int, Int}

  function HcatAtom(args::AbstractExpr...)
    num_rows = args[1].size[1]
    num_cols = 0
    for arg in args
      if arg.size[1] != num_rows
        error("Cannot horizontally stack expressions of varying number of rows")
      end
      num_cols += arg.size[2]
    end
    children = tuple(args...)
    return new(:hcat, hash(children), children, (num_rows, num_cols))
  end
end

function sign(x::HcatAtom)
  return sum(map(sign, x.children))
end

function monotonicity(x::HcatAtom)
  return [Nondecreasing() for c in x.children]
end

function curvature(x::HcatAtom)
  return ConstVexity()
end

function evaluate(x::HcatAtom)
  return hcat([evaluate(c) for c in x.children]...)
end


function conic_form!(x::HcatAtom, unique_conic_forms::UniqueConicForms)
  if !has_conic_form(unique_conic_forms, x)
    # build a list of child conic objectives and constraints
    objectives = ConicObj[]
    for child in x.children
      push!(objectives, conic_form!(child, unique_conic_forms))
    end
    # build a dict from variable ids to sizes
    variable_to_sizes = Dict{UInt64, Int}()
    for objective in objectives
      for id in keys(objective)
        if !(id in keys(variable_to_sizes))
          if id == object_id(:constant)
            variable_to_sizes[id] = 1
          else
            variable_to_sizes[id] = get_vectorized_size(id_to_variables[id])
          end
        end
      end
    end

    # Suppose the child objectives for two children e1 (2 x 1) and e2 (2 x 2) look something like
    #  e1: x => 1 2 3
    #           4 5 6
    #      y => 2 4
    #           7 8
    #  e2: x => 1 1 1
    #           2 2 2
    #           3 3 3
    #           4 4 4
    # The objective of [e1 e2] will look like
    #      x => 1 2 3
    #           4 5 6
    #           1 1 1
    #           2 2 2
    #           3 3 3
    #           4 4 4
    #      y => 2 4
    #           7 8
    #           0 0
    #           0 0
    #           0 0
    #           0 0
    # builds the objective by aggregating a list of coefficients for each variable
    # from each child objective, and then vertically concatenating them
    objective = ConicObj()
    for (id, col_size) in variable_to_sizes
      temp_tuple = Tuple{Value,Value}
      value_list = temp_tuple[]
      for i in 1:length(objectives)
        row_size = get_vectorized_size(x.children[i])
        if haskey(objectives[i], id)
          push!(value_list, objectives[i][id])
        else
          push!(value_list, (spzeros(row_size, col_size),spzeros(row_size, col_size)))
        end
      end
      objective[id] = vcat(value_list...)
    end
    cache_conic_form!(unique_conic_forms, x, objective)
  end
  return get_conic_form(unique_conic_forms, x)
end

hcat(args::AbstractExpr...) = HcatAtom(args...)
hcat(args::AbstractExprOrValue...) = HcatAtom([convert(AbstractExpr, arg) for arg in args]...)
hcat(args::Value...) = Base.cat(2, args...)


# TODO: implement vertical concatenation in a more efficient way
vcat(args::AbstractExpr...) = HcatAtom([arg' for arg in args]...)'
vcat(args::AbstractExprOrValue...) = HcatAtom([convert(AbstractExpr, arg)' for arg in args]...)'
vcat(args::Value...) = Base.cat(1, args...) # Note: this makes general vcat slower for anyone using Convex...
Base.vect{T<:AbstractExpr}(args::T...) = HcatAtom([arg' for arg in args]...)'
Base.vect(args::AbstractExpr...) = HcatAtom([arg' for arg in args]...)'
Base.vect(args::AbstractExprOrValue...) = HcatAtom([convert(AbstractExpr,arg)' for arg in args]...)'
if Base._oldstyle_array_vcat_
    # This is ugly, because the method redefines simple cases like [1,2,3]
    Base.vect(args::Value...) = Base.vcat(args...)
else
    function Base.vect(args::Value...)
        T = Base.promote_typeof(args...)
        return copy!(Array{T}(length(args)), args)
    end
end
