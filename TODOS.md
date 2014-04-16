# Questions

Variable v/s Parameter
Why do variables have linear vexity and what does it mean for a parameter to have a constant vexity
Canonicalization
-what exactly are variables
-we need to create a variable for each node but isnt the AST already doing that

# KV's junk

cd("Desktop/CVX.jl")
include("test/test.jl")

# Next Week's Todos (04/14)

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