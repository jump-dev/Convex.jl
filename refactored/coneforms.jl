function cone_form(x::Constant)
  return ({:Constant, x.value}, ConeConstr[])
end

function cone_form(x::Variable)
  return ({unique_id(x), ones(get_vectorized_size(x))}, ConeConstr[])
end

function cone_form(e::NegateAtom)
  objective, constraints = cone_form(e.children[1])
  for var in keys(objective)
    objective[var] *= -1
  end
  return (objective, constraints)
end

function cone_form(e::AdditionAtom)
  childcones = map(cone_form, e.children)
  objective = Dict()
  constraints = ConeContr[]
  for (childobjective, childconstraints) in childcones
    append!(constraints,childconstraints)
    for var in keys(childobjective)
      objective[var] = has(objective,var) ? childobjective[var] : objective[var] + childobjective[var]
    end
  end
  return (objective, constraints)
end

function cone_form(e::AbsAtom)
  objective, constraints = cone_form(e.children[1])
  t = unique_id(e)
  # t + x >= 0
  constraint1 = copy(objective)
  constraint1[t] = 1
  push!(constraints, ConeConstr(constraint1, :NonNeg))
  # t - x >= 0
  constraint2 = Dict()
  for var in keys(objective)
    constraint2[var] = -objective[var]
  end
  constraint2[t] = 1
  push!(constraints, ConeConstr(constraint2, :NonNeg))
  return ({t=>1}, constraints)  
end
