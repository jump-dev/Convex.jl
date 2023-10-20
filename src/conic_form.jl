# TODO: Comment every single line

# ConicObj represents an affine function of the variables
# it is stored as a dictionary of (key, value) pairs
# keys are unique ids of variables
# values are their coefficients in the affine function
# so for example, {unique_id(x)=>5, unique_id(y)=>6} represents the function 5x + 6y
# we store the affine functions in this form for efficient manipulation of sparse affine functions
struct ConicObj
    mapping::OrderedDict{UInt64,Tuple{Value,Value}}
end
ConicObj() = ConicObj(OrderedDict{UInt64,Tuple{Value,Value}}())
Base.iterate(c::ConicObj, s...) = iterate(c.mapping, s...)
Base.keys(c::ConicObj) = keys(c.mapping)
Base.haskey(c::ConicObj, var::Integer) = haskey(c.mapping, UInt64(var))
Base.getindex(c::ConicObj, var::Integer) = c.mapping[UInt64(var)]
function Base.setindex!(c::ConicObj, val, var::Integer)
    return setindex!(c.mapping, val, UInt64(var))
end
Base.copy(c::ConicObj) = ConicObj(copy(c.mapping))

# helper function to negate conic objectives
# works by changing each (key, val) pair to (key, -val)
function Base.:-(c::ConicObj)
    new_obj = copy(c)
    for var in keys(new_obj)
        x1 = new_obj[var][1] * (-1)
        x2 = new_obj[var][2] * (-1)
        new_obj[var] = (x1, x2)
    end
    return new_obj
end

function Base.:+(c::ConicObj, d::ConicObj)
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
            new_obj[var] = (x1, x2)
        end
    end
    return new_obj
end

function get_row(c::ConicObj, row::Int)
    new_obj = ConicObj()
    for (var, coeff) in c
        x1 = coeff[1][row, :]
        x2 = coeff[2][row, :]
        new_obj[var] = (x1, x2)
    end
    return new_obj
end

function Base.:*(v::Value, c::ConicObj)
    # TODO: this part is time consuming, esp new_obj[var] = v * new_obj[var]...
    new_obj = copy(c)
    for var in keys(new_obj)
        x1 = v * new_obj[var][1]
        x2 = v * new_obj[var][2]
        new_obj[var] = (x1, x2)
    end
    return new_obj
end

function promote_size(c::ConicObj, vectorized_size::Int)
    new_obj = copy(c)
    for var in keys(new_obj)
        x1 = repeat(new_obj[var][1], vectorized_size, 1)
        x2 = repeat(new_obj[var][2], vectorized_size, 1)
        new_obj[var] = (x1, x2)
    end
    return new_obj
end

# A conic constraint is of the form [affine_expr1, affine_expr2, ..., affine_exprk] \in cone
# we represent each affine expressions as a ConicObj
# we represent the cone as a Symbol, like :SOC, :SDP, etc.
# See the details of `get_MOI_set` for the full list of symbols used.
# and we record the sizes of the affine expressions (XXX check...)
# XXX might it be better to represent objs as a single ConicObj rather than an array of them?
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

function Base.empty!(conic_forms::UniqueConicForms)
    empty!(conic_forms.exp_map)
    empty!(conic_forms.constr_map)
    empty!(conic_forms.constr_list)
    empty!(conic_forms.id_to_variables)
    empty!(conic_forms.conic_constr_to_constr)
end

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

function cache_conic_form!(
    conic_forms::UniqueConicForms,
    exp::AbstractExpr,
    new_conic_form::ConicObj,
)
    return conic_forms.exp_map[(exp.head, exp.id_hash)] = new_conic_form
end

function cache_conic_form!(
    conic_forms::UniqueConicForms,
    constr::Constraint,
    new_conic_form::ConicConstr,
)
    conic_forms.constr_map[(constr.head, constr.id_hash)] = 0
    return push!(conic_forms.constr_list, new_conic_form)
end

function cache_conic_form!(
    conic_forms::UniqueConicForms,
    constr::Constraint,
    new_conic_forms::UniqueConstrList,
)
    conic_forms.constr_map[(constr.head, constr.id_hash)] = 0
    return append!(conic_forms.constr_list, new_conic_forms)
end

function add_to_id_to_variables!(
    conic_forms::UniqueConicForms,
    var::AbstractVariable,
)
    return conic_forms.id_to_variables[var.id_hash] = var
end
