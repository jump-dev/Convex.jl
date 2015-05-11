export ConicObj, ConicConstr, UniqueConicForms
export +, -, *, promote_size, get_row
export cache_conic_form!, has_conic_form, get_conic_form

# TODO: Comment every single line
ConicObj = Dict{Uint64, Value}

function -(c::ConicObj)
  new_obj = copy(c)
  for var in keys(new_obj)
    new_obj[var] *= -1
  end
  return new_obj
end

function +(c::ConicObj, d::ConicObj)
  new_obj = copy(c)
  for var in keys(d)
    if !haskey(new_obj, var)
      new_obj[var] = d[var]
    else
      # .+ does not preserve sparsity
      # need to override behavior
      if size(new_obj[var]) == size(d[var])
        new_obj[var] = new_obj[var] + d[var]
      else
        new_obj[var] = broadcast(+, new_obj[var], d[var])
      end
    end
  end
  return new_obj
end

function get_row(c::ConicObj, row::Int)
  new_obj = ConicObj()
  for (var, coeff) in c
    new_obj[var] = coeff[row, :]
  end
  return new_obj
end

function *(v::Value, c::ConicObj)
  # TODO: this part is time consuming, esp new_obj[var] = v * new_obj[var]...
  new_obj = copy(c)
  for var in keys(new_obj)
    new_obj[var] = v * new_obj[var]
  end
  return new_obj
end

function promote_size(c::ConicObj, vectorized_size::Int)
  new_obj = copy(c)
  for var in keys(new_obj)
    new_obj[var] = repmat(new_obj[var], vectorized_size, 1)
  end
  return new_obj
end

type ConicConstr
  objs::Array{ConicObj}
  cone::Symbol
  sizes::Array{Int}
end

UniqueExpMap = Dict{@compat(Tuple{Symbol, Uint64}), ConicObj}
UniqueConstrMap = Dict{@compat(Tuple{Symbol, Uint64}), Int}
UniqueConstrList = Array{ConicConstr}

type UniqueConicForms
  exp_map::UniqueExpMap
  constr_map::UniqueConstrMap
  constr_list::UniqueConstrList
end

UniqueConicForms() = UniqueConicForms(UniqueExpMap(), UniqueConstrMap(), ConicConstr[])

function has_conic_form(conic_forms::UniqueConicForms, exp::AbstractExpr)
  return haskey(conic_forms.exp_map, (exp.head, exp.id_hash))
end

function has_conic_form(conic_forms::UniqueConicForms, constr::Constraint)
  return haskey(conic_forms.constr_map, (constr.head, constr.id_hash))
end

function get_conic_form(conic_forms::UniqueConicForms, exp::AbstractExpr)
  return conic_forms.exp_map[(exp.head, exp.id_hash)]
end

function get_conic_form(conic_forms::UniqueConicForms, constr::Constraint)
  return conic_forms.constr_map[(constr.head, constr.id_hash)]
end

function cache_conic_form!(conic_forms::UniqueConicForms, exp::AbstractExpr, new_conic_form::ConicObj)
  conic_forms.exp_map[(exp.head, exp.id_hash)] = new_conic_form
end

function cache_conic_form!(conic_forms::UniqueConicForms, constr::Constraint, new_conic_form::ConicConstr)
  conic_forms.constr_map[(constr.head, constr.id_hash)] = 0
  push!(conic_forms.constr_list, new_conic_form)
end

function cache_conic_form!(conic_forms::UniqueConicForms, constr::Constraint, new_conic_forms::UniqueConstrList)
  conic_forms.constr_map[(constr.head, constr.id_hash)] = 0
  append!(conic_forms.constr_list, new_conic_forms)
end
