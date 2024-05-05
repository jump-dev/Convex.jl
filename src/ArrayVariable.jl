const GLOBAL_INDEX_COUNTER = Threads.Atomic{UInt64}(0)

struct ScalarVariable{T} <: AbstractVariable
    id_hash::UInt64
    value::Union{T,Nothing}
    function ScalarVariable{T}() where {T}
        # returns the *old* value, so we add 1 to the result
        ret = Threads.atomic_add!(GLOBAL_INDEX_COUNTER, UInt64(1))
        id_hash = Base.checked_add(ret, UInt64(1)) # might as well check...
        return new{T}(id_hash, nothing)
    end
end

# This is like an atom not a variable.
# it is basically a big hvcat of individual scalar variables.
mutable struct ArrayVariable{T, N} <: AbstractExpr
    variables::Array{ScalarVariable{T}, N}
end
function ArrayVariable{T}(size::Tuple) where {T}
    variables = Array{ScalarVariable{T}}(undef, size)
    variables .= ScalarVariable{T}.()
    return ArrayVariable(variables)
end
ArrayVariable{T}(size::Integer...) where {T} = ArrayVariable{T}(size)
ArrayVariable(args...) = ArrayVariable{Float64}(args...)

head(arr::ArrayVariable{T, N}) where {T, N} = "ArrayVariable{$T}($(size(arr)))"
function Base.getproperty(arr::ArrayVariable, field::Symbol)
    if field === :id_hash
        return hash(arr.variables)
    elseif field === :size
        return size(arr)
    elseif field === :value
        return [x.value for x in arr.variables]
    elseif field === :children
        return ()
    else
        return getfield(arr, field)
    end
end

function set_value!(x::ArrayVariable{T}, v::Value) where {T}
    for I in eachindex(v)
        current = x.variables[I]
        x.variables[I] = ScalarVariable{T}(current.id_hash, v[I])
    end
    return nothing
end

vexity(::ArrayVariable) = AffineVexity()
vexity(::ScalarVariable) = AffineVexity()
Base.sign(::ArrayVariable) = NoSign() # real
Base.sign(::ScalarVariable) = NoSign() # real
vartype(::ArrayVariable) = ContVar
vartype(::ScalarVariable) = ContVar
get_constraints(::ScalarVariable) = Constraint[]

Base.size(arr::ArrayVariable) = size(arr.variables)
Base.size(::ScalarVariable) = ()

evaluate(arr::ArrayVariable) = arr.value

const MatrixVariable = ArrayVariable{2}
const VectorVariable = ArrayVariable{1}

export VectorVariable, MatrixVariable, ArrayVariable


function _add_variable(context::Context{T}, a::ScalarVariable) where {T}
    var_index = get!(context.scalar_id_to_moi_index, a.id_hash) do
        return MOI.add_variable(context.model)
    end
    context.id_to_variables[a.id_hash] = a
    return var_index
end

function new_conic_form!(context::Context{T}, arr::ArrayVariable) where {T}
    var_inds = _add_variable.(Ref(context), arr.variables)
    return to_tape(MOI.VectorOfVariables(vec(var_inds)), context)
end

# TODO- indexing, broadcasting, etc
