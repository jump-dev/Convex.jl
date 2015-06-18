# To do

# * consider how to populate dual variables corresponding to eliminated constraints
# * allow users to specify which constraint can't be eliminated
# * populate dual variables for SDPs
# * stop committing crimes against sparse matrices
# * don't give up on eliminating integer / boolean variables
# * eliminate variables and constraints not connected to the objective

# * integers
# * frob norm
# * inv_pos optval

function presolve(c, A, b, constrcones, varcones, vartypes, verbose=false)
	verbose && @show c, full(A), b, constrcones
	m, n = size(A)
	tA = copy(A)
	tb = copy(b)
	tc = copy(c)
	# P*[1; new_primal_optimal_solution] = old_primal_optimal_solution
	P = [eye(n) zeros(n,1)]
	# Don't drop any of the variables present in the varcones, 
	# or integer variables
	keepme = union(append!(Int[vars for (cone,vars) in varcones],
		                  filter(i->vartypes[i]!=:Cont, 1:n)))
	verbose && @show keepme
	eliminatedconstraints = Int[]
	eliminatedvars = Int[]
	
	# we're going to loop over the constraints three times
	# the first time, to detect variables fixed to a given value
	# the second time, to detect equality constraints between just two variables
	# the third time, to delete all-zero rows
	
	# in all other cases we need to eliminate variables in terms of other eliminated variables, 
	# and that sounds dangerously like gaussian elimination, 
	# which we could certainly do much more intelligently.
	# also other cases would cause more fill in

    # note that any variable appearing in a row of A has not previously been eliminated

	# loop 1: detect variables fixed to a given value
	for (cone, indices) in constrcones
		if cone == :Zero			
			for iconstr in indices
				rownnz = find(tA[iconstr,:])
				# row with one entry fixes the value of that variable to be a constant
				if length(rownnz) == 1
					# make sure we don't throw away *all* the variables
					length(eliminatedvars) == n-1 && break
					ivar = rownnz[1]
					ivar in keepme && continue
					P[:,end] += P[:,ivar] * (tb[iconstr] / tA[iconstr,ivar])
					P[:,ivar] = 0
					tb -= tA[:,ivar] * (tb[iconstr] / tA[iconstr,ivar])
					tA[:,ivar] = 0
					tc[ivar] = 0
					push!(eliminatedconstraints, iconstr)
					push!(eliminatedvars, ivar)
				end
			end
		end
	end

	# loop 2: detect equality constraints between just two variables
	for (cone, indices) in constrcones
		if cone == :Zero			
			for iconstr in indices
				rownnz = find(tA[iconstr,:])
				# with two constraints we can always eliminate one variable
				if length(rownnz) == 2
					# we'll eliminate variable with index iv2 and keep iv1
					if rownnz[2] in keepme
						if rownnz[1] in keepme
							continue # don't eliminate either
						end
						# eliminate first variable
						iv1 = rownnz[2]
						iv2 = rownnz[1]
					else
						# eliminate second variable
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
					tb -= Av2*(bc/Acv2)
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
		end
	end
	

	# loop 3: delete all-zero rows
	# we can delete all-zero rows of A in equality and LP constraints
	# and combine them in SOC (as long as it's not the first index) constraints
	# but not in SDP constraints
	# and not in exp constraints
	for (cone, indices) in constrcones
		if cone == :Zero ||	cone == :NonNeg
			for (iconstr,_) in filter(t->length(t[2])==0, map(i->(i,find(tA[i,:])),indices))
				tb[iconstr] >= 0 ? push!(eliminatedconstraints, iconstr) : warn("Infeasible")
			end
		elseif cone == :SOC
			zerorows = filter(t->length(t[2])==0, map(i->(i,find(tA[i,:])),indices[2:end-1]))
			for (iconstr,_) in zerorows[2:end] 
				tb[zerorows[1]] += tb[iconstr]
				push!(eliminatedconstraints, iconstr)
			end
		end
	end

	## Retain only the informative rows and columns of the new problem data
	vars = sort(collect(setdiff(Set(1:n), Set(eliminatedvars))))
	constrs = sort(collect(setdiff(Set(1:m), Set(eliminatedconstraints))))
	verbose && @show vars, constrs

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
		newindices = newconstrindices(indices)
		length(newindices) > 0 && push!(newconstrcones, (cone, newindices))		
	end

	verbose && @show tc[vars], full(tA[constrs, vars]), tb[constrs], newconstrcones, P[:, [vars..., n+1]]
	return tc[vars], tA[constrs, vars], tb[constrs], newconstrcones, vartypes[vars],
		sparse(P[:, push!(vars, n+1)]), constrs
end









