export SOCConstraint, SOCElemConstraint, conic_form

type SOCConstraint <: Constraint
  head::Symbol
  id_hash::Uint64
  children::Tuple

  function SOCConstraint(args::AbstractExpr...)
    children = tuple(args...)
    return new(:soc, hash(children), children)
  end
end

function conic_form(c::SOCConstraint, unique_conic_forms::UniqueConicForms)
  if !has_conic_form(unique_conic_forms, c)
    objectives = Array(ConicObj, length(c.children))
    for iobj=1:length(c.children)
      objectives[iobj] = conic_form(c.children[iobj], unique_conic_forms)
    end
    add_conic_form!(unique_conic_forms, c, ConicConstr(objectives, :SOC, [get_vectorized_size(x) for x in c.children]))
  end
  return get_conic_form(unique_conic_forms, c)
end

type SOCElemConstraint <: Constraint
  head::Symbol
  id_hash::Uint64
  children::Tuple

  function SOCElemConstraint(args::AbstractExpr...)
    children = tuple(args...)
    return new(:soc_elem, hash(children), children)
  end
end

function conic_form(c::SOCElemConstraint, unique_conic_forms::UniqueConicForms)
  if !has_conic_form(unique_conic_forms, c)
    num_constrs = get_vectorized_size(c.children[1])
    num_children = length(c.children)
    conic_constrs = Array(ConicConstr, num_constrs)
    objectives = Array(ConicObj, num_children)
    for iobj = 1:num_children
      objectives[iobj] = conic_form(c.children[iobj], unique_conic_forms)
    end
    for row = 1:num_constrs
      objectives_by_row = [get_row(obj, row) for obj in objectives]
      conic_constrs[row] = ConicConstr(objectives_by_row, :SOC, [1, 1, 1])
    end
    add_conic_form!(unique_conic_forms, c, conic_constrs)
  end
  return get_conic_form(unique_conic_forms, c)
end
