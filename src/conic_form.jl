import Base.+,Base.-,Base.*
export ConicObj, ConicConstr, UniqueConicForms
export +, -, *, promote_size, get_row
export cache_conic_form!, has_conic_form, get_conic_form

# TODO: Comment every single line

# ConicObj represents an affine function of the variables
# it is stored as a disctionary of (key, value) pairs
# keys are unique ids of variables
# values are their coefficients in the affine function
# so for example, {unique_id(x)=>5, unique_id(y)=>6} represents the function 5x + 6y
# we store the affine functions in this form for efficient manipulation of sparse affine functions
ConicObj = DataStructures.OrderedDict{UInt64, Tuple{Value,Value}}

# helper function to negate conic objectives
# works by changing each (key, val) pair to (key, -val)
function -(c::ConicObj)
  new_obj = copy(c)
  for var in keys(new_obj)
    x = new_obj[var][1]*(-1)
    y =  new_obj[var][2]*(-1)
    new_obj[var] = (x,y)
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
      if size(new_obj[var][1]) == size(d[var][1])
        x1 = new_obj[var][1] + d[var][1]
      else
        x1 = broadcast(+, new_obj[var][1], d[var][1])
      end
      if size(new_obj[var][2]) == size(d[var][2])
        x2 = new_obj[var][2] + d[var][2]
      else
        x2 = broadcast(+, new_obj[var][2], d[var][2])
      end
      new_obj[var] = (x1,x2)
    end
  end
  return new_obj
end

function get_row(c::ConicObj, row::Int)
  new_obj = ConicObj()
  for (var, coeff) in c
    x1 = coeff[1][row, :]
    x2 = coeff[2][row, :]
    new_obj[var] = (x1,x2)
  end
  return new_obj
end

function *(v::Value, c::ConicObj)
  # TODO: this part is time consuming, esp new_obj[var] = v * new_obj[var]...
  new_obj = copy(c)
  for var in keys(new_obj)
    x = v * new_obj[var][1]
    y = v * new_obj[var][2]
    new_obj[var] = (x,y)
  end
  return new_obj
end

function promote_size(c::ConicObj, vectorized_size::Int)
  new_obj = copy(c)
  for var in keys(new_obj)
    #for i in 1:2
      x = repmat(new_obj[var][1], vectorized_size, 1)
      y = repmat(new_obj[var][2], vectorized_size, 1)
      new_obj[var] = (x,y)
    #end
  end
  return new_obj
end

# A conic constraint is of the form [affine_expr1, affine_expr2, ..., affine_exprk] \in cone
# we represent each affine expressions as a ConicObj
# we represent the cone as a Symbol (defined in MathProgBase), like :SOC, :LP, etc
# and we record the sizes of the affine expressions (XXX check...)
# XXX might it be better to represent objs as a single ConicObj rather than an array of them?
type ConicConstr
  objs::Array{ConicObj}
  cone::Symbol
  sizes::Array{Int}
end

# in conic form, every expression e is represented by a ConicObj together with a collection of ConicConstrs
# for each expression e, UniqueExpMap maps (e.head, unique_id(e)) to that expression's ConicObj 
UniqueExpMap = DataStructures.OrderedDict{Tuple{Symbol, UInt64}, ConicObj}
# for each expression e, UniqueExpMap maps (e.head, unique_id(e)) to the index of expression's ConicConstr in UniqueConstrList
UniqueConstrMap = DataStructures.OrderedDict{Tuple{Symbol, UInt64}, Int}
# records each ConicConstr created
UniqueConstrList = Array{ConicConstr}

# UniqueConicForms caches all the conic forms of expressions we've parsed so far
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
