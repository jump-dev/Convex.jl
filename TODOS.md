point of zero sign

What happens if you take something which is convex and multiply it with something with sign any

max and stuff

do we worry about parameters now

getting dual variables?

The constraints list that we pass in, each constraint has some uid. So we just need to find the duals associated with those constraints (not the constraints created by all the other variables we created in canonicalization)

We should focus on LPs, but the only reason we might want to focus on SOCPs is because there are no conic solvers in Julia

1. delay computation of dual variables
2. get a a wrapper to scs and figure out a canon form for cone problems


TODOS
1. wrapper for SCS
2. we can look at GLPK. One time madeleine wanted a corner point so she used a simplex solver rather than an interior point solver. something something? what? something that needs to be done but not urgent?


# Next Week's Todos (04/14)

- Next week, we want to finish standard form LPs completely.
- Implement canonicalization: Hopefully there won't be big efficiency concerns here, it's mostly just up to how the language is implemented.
- Add tests
- Read Performance Tips from Julia docs

# This Week's Updates

- Baby ECOS problem

---

# Big Picture

For our EE 364B project, we aim to implement an LP solver for CVX.jl. There are some basic steps (required ones in bold)

1. **Abstract syntax tree**
2. (Optional) Reductions to make it a simpler abstract syntax tree
3. **Canonicalization**: Introduce new variables
4. _Relaxation_: (Not necessary for LPs) For variables that come from SOCP or SDP constraints, relax to a convex program
5. **Stuff the matrix _symbolically_**: How to differentiate CVX.jl from existing things like CVXPY
6. Turn the block matrix into a literal / **numerical matrix**
7. **Call a solver**: ECOS or SCS

# Notes

What are the bottlenecks in existing solvers like CVXPY?

- Declaring new variables: Julia is less object oriented so it might be better. Python has slow loops, Julia has fast loops.
- Copying all the data
