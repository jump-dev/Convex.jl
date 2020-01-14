#############################################################################
# variable.jl
# Defines Variable, which is a subtype of AbstractExpr
#############################################################################

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


function vexity(x::Variable)
    return x.vexity
end

function evaluate(x::Variable)
    return x.value === nothing ? error("Value of the variable is yet to be calculated") : output(x.value)
end

function sign(x::Variable)
    return x.sign
end


function real_conic_form(x::Variable)
    vec_size = length(x)
    return sparse(1.0I, vec_size, vec_size)
end

function imag_conic_form(x::Variable)
    vec_size = length(x)
    if x.sign == ComplexSign()
        return im*sparse(1.0I, vec_size, vec_size)
    else
        return spzeros(vec_size, vec_size)
    end
end

# fix variables to hold them at their current value, and free them afterwards
function fix!(x::Variable)
    x.value === nothing && error("This variable has no value yet; cannot fix value to nothing!")
    x.vexity = ConstVexity()
    x
end
function fix!(x::Variable, v::AbstractArray)
    size(x) == size(v) || throw(DimensionMismatch("Variable and value sizes do not match!"))
    x.value = sign(x) == ComplexSign() ? convert(Array{ComplexF64}, v) : convert(Array{Float64}, v)
    fix!(x)
end

function fix!(x::Variable, v::AbstractVector)
    size(x, 2) == 1 || throw(DimensionMismatch("Cannot set value of a variable of size $(size(x)) to a vector"))
    size(x, 1) == length(v) || throw(DimensionMismatch("Variable and value sizes do not match!"))
    x.value = sign(x) == ComplexSign() ? convert(Array{ComplexF64}, v) : convert(Array{Float64}, v)
    fix!(x)
end

function fix!(x::Variable, v::Number)
    size(x) == (1,1) || throw(DimensionMismatch("Variable and value sizes do not match!"))
    x.value = sign(x) == ComplexSign() ? convert(ComplexF64, v) : convert(Float64, v)
    fix!(x)
end

function free!(x::Variable)
    x.vexity = AffineVexity()
    x
end
