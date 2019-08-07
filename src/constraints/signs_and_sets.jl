export conic_form!

conic_form!(s::Positive, x::AbstractVariable, unique_conic_forms) = conic_form!(x>=0, unique_conic_forms)
conic_form!(s::Negative, x::AbstractVariable, unique_conic_forms) = conic_form!(x<=0, unique_conic_forms)