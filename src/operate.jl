# Copyright (c) 2014: Madeleine Udell and contributors
#
# Use of this source code is governed by a BSD-style license that can be found
# in the LICENSE file or at https://opensource.org/license/bsd-2-clause

complex_promote(tape::ComplexTape) = tape
complex_promote(tape::SparseTape) = tape
complex_promote(v::ComplexStructOfVec) = v
complex_promote(x::Real) = complex(x)
complex_promote(x::Complex) = x

function complex_promote(v::AbstractVector{T}) where {T}
    return ComplexStructOfVec(v, spzeros(T, length(v)))
end

# Here we run `complex_promote` and dispatch to either `real_operate`
# or `complex_operate` depending on `sign`. We also check the types
# at this stage, instead of letting them fall to MethodErrors later.
# At this stage, it should be impossible for any user input to lead to an
# error here, since we have converted the arguments earlier in the pipeline.
# However new atoms and new invocations of `operate` could lead to issues.
function operate(op::F, ::Type{T}, sign::Sign, args...) where {F,T}
    if iscomplex(sign)
        if op === add_operation
            @assert length(args) == 2
            x, y = args
            if !(
                typeof(x) in (
                    T,
                    Complex{T},
                    SPARSE_MATRIX{T},
                    SPARSE_MATRIX{Complex{T}},
                    Matrix{T},
                    Matrix{Complex{T}},
                )
            )
                error(
                    "Convex.jl internal error: unexpected type $(typeof(x)) for first argument of operation $op of with sign=$sign",
                )
            end
            if !(
                typeof(y) in (
                    SparseTape{T},
                    SPARSE_VECTOR{T},
                    ComplexTape{T},
                    ComplexStructOfVec{T},
                )
            )
                error(
                    "Convex.jl internal error: unexpected type $(typeof(y)) for second argument of operation $op of with sign=$sign",
                )
            end
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
                    error(
                        "Convex.jl internal error: unexpected type $(typeof(arg))",
                    )
                end
            end
        end
        return complex_operate(op, T, args...)
    else
        if op in (real, imag, add_operation)
            if op in (real, imag)
                x = only(args)
                x = complex_promote(x)

                if !(
                    typeof(x) == ComplexTape{T} ||
                    typeof(x) == ComplexStructOfVec{T} ||
                    typeof(x) == SparseTape{T}
                )
                    error(
                        "Convex.jl internal error: unexpected type $(typeof(x))",
                    )
                end
            end

            if op === add_operation
                @assert length(args) == 2

                x, y = args
                if !(typeof(x) in (T, SPARSE_MATRIX{T}, Matrix{T}))
                    error(
                        "Convex.jl internal error: unexpected type $(typeof(x)) for first argument of operation $op of with sign=$sign",
                    )
                end
                if !(typeof(y) in (SparseTape{T}, SPARSE_VECTOR{T}))
                    error(
                        "Convex.jl internal error: unexpected type $(typeof(y)) for second argument of operation $op of with sign=$sign",
                    )
                end
            end
        else
            # Everything should be either a
            # SparseTape{T} or SPARSE_VECTOR{T}
            for arg in args
                if !(
                    typeof(arg) == SparseTape{T} ||
                    typeof(arg) == SPARSE_VECTOR{T}
                )
                    error("Internal error: unexpected type $(typeof(arg))")
                end
            end
        end
        return real_operate(op, T, args...)
    end
end
