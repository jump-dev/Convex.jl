function conic_form!(x::Variable, unique_conic_forms::UniqueConicForms=UniqueConicForms())
    if !has_conic_form(unique_conic_forms, x)
        cache_conic_form!(unique_conic_forms, x)
        if vexity(x) == ConstVexity()
            # do exactly what we would for a constant
            objective = ConicObj()
            objective[objectid(:constant)] = (vec([real(x.value);]),vec([imag(x.value);]))
            cache_conic_form!(unique_conic_forms, x, objective)
        else
            objective = ConicObj()
            vec_size = length(x)

            objective[x.id_hash] = (real_conic_form(x), imag_conic_form(x))
            objective[objectid(:constant)] = (spzeros(vec_size, 1), spzeros(vec_size, 1))
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