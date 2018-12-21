#############################################################################
# variable.jl
# Defines Variable, which is a subtype of AbstractExpr
#############################################################################

export Variable, Semidefinite, ComplexVariable, HermitianSemidefinite
export vexity, evaluate, sign, conic_form!, fix!, free!

mutable struct Variable <: AbstractExpr
    head::Symbol
    id_hash::UInt64
    value::ValueOrNothing
    size::Tuple{Int, Int}
    vexity::Vexity
    sign::Sign
    sets::Array{Symbol,1}


    function Variable(size::Tuple{Int, Int}, sign::Sign=NoSign(), sets::Symbol...)
        this = new(:variable, 0, nothing, size, AffineVexity(), sign, Symbol[sets...])
        this.id_hash = objectid(this)
        id_to_variables[this.id_hash] = this
        return this
    end

    Variable(m::Int, n::Int, sign::Sign=NoSign(), sets::Symbol...) = Variable((m,n), sign, sets...)
    Variable(sign::Sign, sets::Symbol...) = Variable((1, 1), sign, sets...)
    Variable(sets::Symbol...) = Variable((1, 1), NoSign(), sets...)
    Variable(size::Tuple{Int, Int}, sets::Symbol...) = Variable(size, NoSign(), sets...)
    Variable(size::Int, sign::Sign=NoSign(), sets::Symbol...) = Variable((size, 1), sign, sets...)
    Variable(size::Int, sets::Symbol...) = Variable((size, 1), sets...)
end

Semidefinite(m::Integer) = Variable((m, m), :Semidefinite)
function Semidefinite(m::Integer, n::Integer)
    if m == n
        return Variable((m, m), :Semidefinite)
    else
        error("Semidefinite matrices must be square")
    end
end

ComplexVariable(m::Int, n::Int, sets::Symbol...) = Variable((m, n), ComplexSign(), sets...)
ComplexVariable(sets::Symbol...) = Variable((1, 1), ComplexSign(), sets...)
ComplexVariable(size::Tuple{Int, Int}, sets::Symbol...) = Variable(size, ComplexSign(), sets...)
ComplexVariable(size::Int, sets::Symbol...) = Variable((size, 1), ComplexSign(), sets...)
HermitianSemidefinite(m::Integer) = ComplexVariable((m, m), :Semidefinite)
function HermitianSemidefinite(m::Integer, n::Integer)
    if m == n
        return ComplexVariable((m, m), :Semidefinite)
    else
        error("HermitianSemidefinite matrices must be square")
    end
end

# global map from unique variable ids to variables.
# the expression tree will only utilize variable ids during construction
# full information of the variables will be needed during stuffing
# and after solving to populate the variables with values
id_to_variables = Dict{UInt64, Variable}()

function vexity(x::Variable)
    return x.vexity
end

function evaluate(x::Variable)
    return x.value === nothing ? error("Value of the variable is yet to be calculated") : x.value
end

function sign(x::Variable)
    return x.sign
end


function real_conic_form(x::Variable)
    return Eye{Float64}(length(x))
end

function imag_conic_form(x::Variable)
    vec_size = length(x)
    if x.sign == ComplexSign()
        return im * Eye{Float64}(vec_size)
    else
        return Zeros{Float64}(vec_size, vec_size)
    end
end

function conic_form!(x::Variable, unique_conic_forms::UniqueConicForms=UniqueConicForms())
    if !has_conic_form(unique_conic_forms, x)
        if :fixed in x.sets
            # do exactly what we would for a constant
            objective = ConicObj()
            objective[objectid(:constant)] = (vec([real(x.value);]),vec([imag(x.value);]))
            cache_conic_form!(unique_conic_forms, x, objective)
        else
            objective = ConicObj()
            vec_size = length(x)

            objective[x.id_hash] = (real_conic_form(x), imag_conic_form(x))
            objective[objectid(:constant)] = (Zeros{Float64}(vec_size, 1), Zeros{Float64}(vec_size, 1))
            # placeholder values in unique constraints prevent infinite recursion depth
            cache_conic_form!(unique_conic_forms, x, objective)
            if !(x.sign == NoSign() || x.sign == ComplexSign())
                conic_form!(x.sign, x, unique_conic_forms)
            end
            for set in x.sets
                conic_form!(set, x, unique_conic_forms)
            end
        end
    end
    return get_conic_form(unique_conic_forms, x)
end

# fix variables to hold them at their current value, and free them afterwards
function fix!(x::Variable)
    x.value === nothing && error("This variable has no value yet; cannot fix value to nothing!")
    push!(x.sets, :fixed)
    x.vexity = ConstVexity()
    x
end
function fix!(x::Variable, v::AbstractArray)
    size(x) == size(v) || throw(DimensionMismatch("Variable and value sizes do not match!"))
    x.value = Array{Float64}(undef, size(x))
    x.value[:,:] = v
    fix!(x)
end
fix!(x::Variable, v::Number) = fix!(x, fill(v, (1, 1)))

function free!(x::Variable)
    # TODO this won't work if :fixed appears other than at the end of x.sets
    x.sets[end] == :fixed && pop!(x.sets)
    x.vexity = AffineVexity()
    x
end
