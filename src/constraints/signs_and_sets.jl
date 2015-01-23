export conic_form!

conic_form!(s::Positive, x::Variable, unique_conic_forms) = conic_form!(x>=0, unique_conic_forms)
conic_form!(s::Negative, x::Variable, unique_conic_forms) = conic_form!(x<=0, unique_conic_forms)

function conic_form!(set::Symbol, x::Variable, unique_conic_forms)
	if set==:Semidefinite
		conic_form!(SDPConstraint(x), unique_conic_forms)
	end
end
