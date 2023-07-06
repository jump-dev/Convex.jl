# TODO: Delete

# ConicObj represents an affine function of the variables
# it is stored as a dictionary of (key, value) pairs
# keys are unique ids of variables
# values are their coefficients in the affine function
# so for example, {unique_id(x)=>5, unique_id(y)=>6} represents the function 5x + 6y
# we store the affine functions in this form for efficient manipulation of sparse affine functions
struct ConicObj
    mapping::OrderedDict{UInt64,Tuple{Value,Value}}
end
struct ConicConstr
    objs::Vector{ConicObj}
    cone::Symbol
    sizes::Vector{Int}
end

# in conic form, every expression e is represented by a ConicObj together with a collection of ConicConstrs
# for each expression e, UniqueExpMap maps (e.head, unique_id(e)) to that expression's ConicObj
const UniqueExpMap = OrderedDict{Tuple{Symbol,UInt64},ConicObj}
# for each expression e, UniqueExpMap maps (e.head, unique_id(e)) to the index of expression's ConicConstr in UniqueConstrList
const UniqueConstrMap = OrderedDict{Tuple{Symbol,UInt64},Int}
# records each ConicConstr created
const UniqueConstrList = Vector{ConicConstr}
# map variables' hash to the variable itself
const IdToVariables = OrderedDict{UInt64,AbstractVariable}
const ConicConstrToConstr = OrderedDict{ConicConstr,Constraint}
# UniqueConicForms caches all the conic forms of expressions we've parsed so far
struct UniqueConicForms
    exp_map::UniqueExpMap
    constr_map::UniqueConstrMap
    constr_list::UniqueConstrList
    id_to_variables::IdToVariables
    conic_constr_to_constr::ConicConstrToConstr
end

function UniqueConicForms()
    return UniqueConicForms(
        UniqueExpMap(),
        UniqueConstrMap(),
        ConicConstr[],
        IdToVariables(),
        ConicConstrToConstr(),
    )
end
