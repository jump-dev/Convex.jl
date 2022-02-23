function conic_form!(s::Positive, x::AbstractVariable, unique_conic_forms)
    return conic_form!(x >= 0, unique_conic_forms)
end
function conic_form!(s::Negative, x::AbstractVariable, unique_conic_forms)
    return conic_form!(x <= 0, unique_conic_forms)
end
