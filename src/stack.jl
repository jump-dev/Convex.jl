import Base.vcat, Base.hcat
export vcat, hcat, VcatAtom, HcatAtom
export sign, curvature, monotonicity, evaluate, conic_form

type HcatAtom <: AbstractExpr
  head::Symbol
  children_hash::Uint64
  children::Tuple
  size::(Int64, Int64)

  function HcatAtom(args::AbstractExpr...)
    num_rows = args[1].size[1]
    num_cols = 0
    for arg in args
      if arg.size[1] != num_rows
        error("Cannot horizonatally stack expressions of varying number of rows")
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

# TODO: evaluate


function conic_form(x::HcatAtom, unique_constr)
  if !((x.head, x.children_hash) in unique_constr)
    # build a list of child conic objectives and constraints
    constraints = ConicConstr[]
    objectives = ConicObj[]
    for child in x.children
      child_obj, child_constrs = conic_form(child)
      push!(objectives, child_obj)
      append!(constraints, child_constrs)
    end

    # build a dict from variable ids to coefficient values
    variable_to_values = Dict{Uint64, Array{Value}}()
    for objective in objectives
      for (id, value) in objective
        if !(id in variable_to_values)
          variable_to_values[id] = Value[]
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
    num_rows = x.size[1]
    num_cols = x.size[2]
    for (id, value_list) in variable_to_values
      var_size = get_vectorized_size(id_to_variables[id])
      for i in 1:length(objectives)
        if id in objectives[i]
          push!(value_list, objectives[i][id])
        else
          push!(value_list, spzeros(num_rows * num_cols, var_size))
        end
      end
      objective[id] = vcat(value_list...)
    end
    unique_constr[(x.head, x.children_hash)] = (objective, constraints)
  end
  return unique_constr[(x.head, x.children_hash)]
end

hcat(args::AbstractExpr...) = HcatAtom(args...)



type VcatAtom <: AbstractExpr
  head::Symbol
  children_hash::Uint64
  chilren::Tuple
  size::(Int64, Int64)

  function VcatAtom(args::AbstractExpr...)
    num_rows = 0
    num_cols = args[1].size[2]
    for arg in args
      if arg.size[2] != num_cols
        error("Cannot vertically stack expressions of varying number of columns")
      end
      num_rows += args.size[1]
    end
    children = tuple(args...)
    return new(:vcat, hash(children), children, (num_rows, num_cols))
  end
end

function sign(x::VcatAtom)
  return sum(map(sign, x.children))
end

function monotonicity(x::VcatAtom)
  return [Nondecreasing() for c in x.children]
end

function curvature(x::VcatAtom)
  return ConstVexity()
end


# TODO evaluate

function conic_form(x::VcatAtom, unique_constr)
  # build a list of child conic objectives and constraints
  constraints = ConicConstr[]
  objectives = ConicObj[]
  for child in x.children
    child_obj, child_constrs = conic_form(child)
    push!(objectives, child_obj)
    push!(constraints, child_constrs)
  end

  # build a dict from variable ids to coefficient values
  variable_to_values = Dict{Uint64, Value[]}()
  for objective in objectives
    for (id, value) in objective
      if !(id in variable_to_values)
        variable_to_values[id] = Value[]
      end
    end
  end

  # Suppose the child objectives for two children e1 (1 x 2) and e2 (2 x 2) look something like
  #  e1: x => 1 2 3
  #           4 5 6
  #      y => 2 4
  #           7 8
  #  e2: x => 1 1 1
  #           2 2 2
  #           3 3 3
  #           4 4 4
  # The objective of [e1, e2] will look like
  #      x => 1 2 3
  #           1 1 1
  #           2 2 2
  #           4 5 6
  #           3 3 3
  #           4 4 4
  #      y => 2 4
  #           0 0
  #           0 0
  #           7 8
  #           0 0
  #           0 0

  # TODO how to best implement the above?

end
