# To do

# * consider how to populate dual variables corresponding to eliminated constraints
# * allow users to specify which constraint can't be eliminated
# * populate dual variables for SDPs
# * stop committing crimes against sparse matrices

function presolve(c, A, b, constrcones, varcones)
	@show c, full(A), b, constrcones
	m, n = size(A)
	tA = copy(A)
	tb = copy(b)
	tc = copy(c)
	# P*[1; new_primal_optimal_solution] = old_primal_optimal_solution
	P = [eye(n) zeros(n,1)]
	# Don't drop any of the variables present in the varcones
	keepme = union([vars for (cone,vars) in varcones])
	eliminatedconstraints = Int[]
	eliminatedvars = Int[]
	for (cone, indices) in constrcones
		if cone == :Zero			
			for iconstr in indices
				rownnz = find(tA[iconstr,:])
				# note that any variable appearing in a row of A has not previously been eliminated
				# all zero row can be eliminated
				if length(rownnz) == 0
					tb[iconstr] == 0 ? push!(eliminatedconstraints, iconstr) : error("Infeasible")
				# row with one entry fixes the value of that variable to be a constant
				elseif length(rownnz) == 1
					# make sure we don't throw away *all* the variables
					length(eliminatedvars) == n-1 && break
					ivar = rownnz[1]
					P[:,end] += P[:,ivar] * (tb[iconstr] / tA[iconstr,ivar])
					P[:,ivar] = 0
					tb += tA[:,ivar] * (tb[iconstr] / tA[iconstr,ivar])
					tA[:,ivar] = 0
					tc[ivar] = 0
					push!(eliminatedconstraints, iconstr)
					push!(eliminatedvars, ivar)
				# with two constraints we can always eliminate one variable
				elseif length(rownnz) == 2
					if !(rownnz[2] in keepme)
						# eliminate second variable
						iv1 = rownnz[1]
						iv2 = rownnz[2]
					else
						# eliminate first variable
						iv1 = rownnz[1]
						iv2 = rownnz[2]
					end
					bc = tb[iconstr]
					Acv1 = tA[iconstr, iv1]
					Acv2 = tA[iconstr, iv2]
					Av2 = tA[:,iv2]
					cv1 = tc[iv1]
					cv2 = tc[iv2]
					# Acv1*v1 + Acv2*v2 = bc => v2 = bc/Acv2 - Acv1/Acv2*v1
					P[iv2, end] = bc / Acv2
					P[iv2, iv2] = 0
					P[iv2, iv1] = - Acv1/Acv2
					# Av2*v2 = Av2*(bc/Acv2 - Acv1/Acv2*v1)
					tb += Av2*(bc/Acv2)
					tA[:,iv1] -= Av2*(Acv1/Acv2)
					tA[:,iv2] = 0
					# cv2*v2 = cv2*(bc/Acv2 - Acv1/Acv2*v1), but ignore the constant term
					tc[iv1] -= cv2*Acv1/Acv2
					tc[iv2] = 0
					push!(eliminatedconstraints, iconstr)
					push!(eliminatedvars, iv2)
					push!(eliminatedconstraints, iconstr)
				end
			end
			# in all other cases we need to eliminate variables in terms of other eliminated variables, 
			# and that sounds dangerously like gaussian elimination, 
			# which we could certainly do much more intelligently.
			# also other cases would cause more fill in
		end
		# btw we can also delete all-zero rows in LP and SOC constraints.
		# not in SDP constraints.
		# and this process might introduce new all-zero rows.
		# so we could check for them in a second loop at the end if that seemed important.
	end

	## Retain only the informative rows and columns of the new problem data
	vars = sort(collect(union(keepme, setdiff(Set(1:n), Set(eliminatedvars)))))
	constrs = sort(collect(setdiff(Set(1:m), Set(eliminatedconstraints))))
	@show vars, constrs

	# map constraint indices onto new constraint indices in constraintcones
	constrmap = zeros(Int, m)
	iconstr = 1
	for i=1:m
		if !(i in eliminatedconstraints)
			constrmap[i] = iconstr
			iconstr += 1
		end
	end
	newconstrindices(indices) = filter(x->x>0, constrmap[indices])
	newconstrcones = (Symbol, Any)[]
	for (cone,indices) in constrcones
		push!(newconstrcones, (cone, newconstrindices(indices)) )		
	end

	@show tc[vars], full(tA[constrs, vars]), tb[constrs], newconstrcones, 
		P[:, [vars..., n+1]]
	return tc[vars], tA[constrs, vars], tb[constrs], newconstrcones, 
		sparse(P[:, push!(vars, n+1)]), constrs
end









