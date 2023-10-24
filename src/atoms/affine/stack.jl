import Base.vcat, Base.hcat, Base.hvcat
struct HcatAtom <: AbstractExpr
    head::Symbol
    id_hash::UInt64
    children::Tuple
    size::Tuple{Int,Int}

    function HcatAtom(args)
        if !all(x -> x isa AbstractExpr, args)
            args = map(arg -> convert(AbstractExpr, arg), args)
        end
        num_rows = args[1].size[1]
        num_cols = 0
        for arg in args
            if arg.size[1] != num_rows
                error(
                    "Cannot horizontally stack expressions of varying number of rows",
                )
            end
            num_cols += arg.size[2]
        end
        children = tuple(args...)
        return new(:hcat, hash(children), children, (num_rows, num_cols))
    end
end
# Allow splatting
HcatAtom(args...) = HcatAtom(args)

function sign(x::HcatAtom)
    return sum(map(sign, x.children))
end

function monotonicity(x::HcatAtom)
    return [Nondecreasing() for c in x.children]
end

function curvature(x::HcatAtom)
    return ConstVexity()
end

function evaluate(x::HcatAtom)
    return hcat(map(evaluate, x.children)...)
end

function conic_form!(x::HcatAtom, unique_conic_forms::UniqueConicForms)
    if !has_conic_form(unique_conic_forms, x)
        # build a list of child conic objectives and constraints
        objectives = ConicObj[]
        for child in x.children
            push!(objectives, conic_form!(child, unique_conic_forms))
        end
        # build a dict from variable ids to sizes
        variable_to_sizes = OrderedDict{UInt64,Int}()
        for objective in objectives
            for id in keys(objective)
                if !(id in keys(variable_to_sizes))
                    if id == objectid(:constant)
                        variable_to_sizes[id] = 1
                    else
                        variable_to_sizes[id] =
                            length(unique_conic_forms.id_to_variables[id])
                    end
                end
            end
        end

        # Suppose the child objectives for two children e1 (2 x 1) and e2 (2 x 2) look something like
        #  e1: x => 1 2 3
        #           4 5 6
        #      y => 2 4
        #           7 8
        #  e2: x => 1 1 1
        #           2 2 2
        #           3 3 3
        #           4 4 4
        # The objective of [e1 e2] will look like
        #            x => 1 2 3
        #                 4 5 6
        #                 1 1 1
        #                 2 2 2
        #                 3 3 3
        #                 4 4 4
        #            y => 2 4
        #                 7 8
        #                 0 0
        #                 0 0
        #                 0 0
        #                 0 0
        # builds the objective by aggregating a list of coefficients for each variable
        # from each child objective, and then vertically concatenating them
        objective = ConicObj()
        for (id, col_size) in variable_to_sizes
            #temp_tuple = Tuple{Value,Value}
            x1_value_list = Value[]
            x2_value_list = Value[]

            for i in 1:length(objectives)
                row_size = length(x.children[i])
                if haskey(objectives[i], id)
                    push!(x1_value_list, objectives[i][id][1])
                    push!(x2_value_list, objectives[i][id][2])
                else
                    push!(x1_value_list, spzeros(row_size, col_size))
                    push!(x2_value_list, spzeros(row_size, col_size))
                end
            end
            x1 = vcat(x1_value_list...)
            x2 = vcat(x2_value_list...)
            objective[id] = (x1, x2)
        end
        cache_conic_form!(unique_conic_forms, x, objective)
    end
    return get_conic_form(unique_conic_forms, x)
end

# Avoid piracy: if any 1 argument is an AbstractExpr, out of up to 5 total arguments
# then dispatch to our private method (e.g. `HcatAtom` for `hcat`).
# Uses the `eval` strategy "At least one argument of this type" from
# from https://www.oxinabox.net/2023/06/17/resident-eval.html
# How to choose how many arguments to support?
# - Up to 5 arguments means 57 new methods for each function
# - Up to 6 arguments means 120 new methods for each function
# - Up to 7 arguments means 247 new methods for each function
const N_METHODS = @load_preference("num_cat_methods", 5)

# Let us stick to 5 unless we really need more.
for (outer, inner) in [(:hcat, :HcatAtom), (:vcat, :_vcat)]
    for len in 1:N_METHODS  # generate all combinations up to length 5
        for mask in Iterators.product(ntuple(_ -> (true, false), len)...)
            any(mask) || continue  # Don't do this if no argument would be a Foo
            arg_names = Symbol[]
            sig = Expr[]
            for (ii, is_foo) in enumerate(mask)
                arg_name = Symbol(:x, ii)
                push!(arg_names, arg_name)
                push!(sig, :($arg_name::$(is_foo ? :AbstractExpr : :Value)))
            end
            body = quote
                $inner($(arg_names...))
            end
            eval(Expr(:function, Expr(:call, outer, sig...), body))
        end
    end
end

# `hvcat` is special since the first argument is different
for len in 1:N_METHODS  # generate all combinations up to length 5
    for mask in Iterators.product(ntuple(_ -> (true, false), len)...)
        any(mask) || continue  # Don't do this if no argument would be a Foo
        arg_names = Symbol[:rows]
        sig = Expr[:(rows::Tuple{Vararg{Int}})]
        for (ii, is_foo) in enumerate(mask)
            arg_name = Symbol(:x, ii)
            push!(arg_names, arg_name)
            push!(sig, :($arg_name::$(is_foo ? :AbstractExpr : :Value)))
        end
        body = quote
            _hvcat($(arg_names...))
        end
        eval(Expr(:function, Expr(:call, :hvcat, sig...), body))
    end
end

# TODO: implement vertical concatenation in a more efficient way
_vcat(args::AbstractExpr...) = transpose(HcatAtom(map(transpose, args)))

function _vcat(args::AbstractExprOrValue...)
    return transpose(
        HcatAtom(map(arg -> transpose(convert(AbstractExpr, arg)), args)),
    )
end

function _hvcat(rows::Tuple{Vararg{Int}}, args::AbstractExprOrValue...)
    nbr = length(rows)
    rs = Vector{Any}(undef, nbr)
    a = 1
    for i in 1:nbr
        # Avoid splatting to not hit method limits
        rs[i] = HcatAtom(args[a:a-1+rows[i]])
        a += rows[i]
    end
    return vcat(rs...)
end
