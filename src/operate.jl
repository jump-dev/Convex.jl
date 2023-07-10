complex_promote(tape::ComplexTape) = tape
complex_promote(tape::SparseTape) = tape
complex_promote(v::ComplexStructOfVec) = v

function complex_promote(v::Vector)
    return ComplexStructOfVec(v, zero(v))
end

mapping = Dict()
function operate(op::F, ::Type{T}, sign::Sign, args...) where {F,T}
    v = get!(mapping, (op, iscomplex(sign)), Any[])
    push!(v, map(typeof, args))

    if iscomplex(sign)
        if op === add_operation

        else
            # Everything should be either be
            # ComplexTape{T}, SparseTape{T}, or ComplexStructOfVec{T}
            args = map(complex_promote, args)
            for arg in args
                if !(
                    typeof(arg) == ComplexTape{T} ||
                    typeof(arg) == ComplexStructOfVec{T} ||
                    typeof(arg) == SparseTape{T}
                )
                    error("Internal error: unexpected type $(typeof(arg))")
                end
            end
        end
        return complex_operate(op, T, args...)
    else
        if op in (real, imag, add_operation)
            # TODO- add some kind of checks
            # @assert length(args) == 2
            # x, y = args
            # if !(
            #     x isa Union{T,SparseMatrixCSC{T}} && y isa SparseTape{T} ||
            #     y isa Union{T,SparseMatrixCSC{T}} && x isa SparseTape{T}
            # )
            #     error(
            #         "Internal error: unexpected type pair $(typeof(x)) and $(typeof(y))",
            #     )
            # end
        else
            # Everything should be either a
            # SparseTape{T} or Vector{T}
            # Except for `real` and `imag` and `add_operation`...
            for arg in args
                if !(typeof(arg) == SparseTape{T} || typeof(arg) == Vector{T})
                    error("Internal error: unexpected type $(typeof(arg))")
                end
            end
        end
        return real_operate(op, T, args...)
    end
end
