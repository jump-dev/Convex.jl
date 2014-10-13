export ConicObj, ConicConstr
export +, -, *, promote_size, safe_copy

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

function *(v::Value, c::ConicObj)
  new_obj = copy(c)
  for var in keys(new_obj)
    new_obj[var] = v * new_obj[var]
  end
  return new_obj
end

function promote_size(c::ConicObj, vectorized_size::Int64)
  new_obj = copy(c)
  for var in keys(new_obj)
    new_obj[var] = repmat(new_obj[var], vectorized_size, 1)
  end
  return new_obj
end

type ConicConstr
  objs::Array{ConicObj}
  cone::Symbol
  sizes::Array{Int64}
end

function safe_copy(c::(ConicObj, Array{ConicConstr}))
  return deepcopy(c[1]), deepcopy(c[2])
end
