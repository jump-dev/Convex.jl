export ConicObj, ConicConstr
export +

# TODO: @david, just make ConicObj = Dict{Uint64, Value}
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
      # .+ does not behave properly for sparse matrices
      # need to override behavior
      new_obj[var] += d[var]
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


type ConicConstr
  objs::Array{ConicObj}
  cone::Symbol
  sizes::Array{Int64}
end

