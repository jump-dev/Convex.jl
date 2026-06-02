# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Convex.jl is a Julia package for Disciplined Convex Programming (DCP). Users describe optimization problems in natural mathematical syntax, and the package reformulates them into conic programs solved via MathOptInterface (MOI)-compatible solvers (SCS, Mosek, ECOS, Clarabel, etc.).

## Common Commands

### Testing
```bash
# Run full test suite
julia --project=test -e "using Pkg; Pkg.test()"

# Run a specific test module (atoms, utilities, constraints, problem_depot, MOI_wrapper)
julia --project=test -e "include(\"test/test_atoms.jl\"); TestAtoms.runtests()"
julia --project=test -e "include(\"test/test_utilities.jl\"); TestUtilities.runtests()"
julia --project=test -e "include(\"test/test_constraints.jl\"); TestConstraints.runtests()"
```

### Formatting
Code must pass JuliaFormatter checks (margin=80, `always_for_in=true`, `always_use_return=true`, `remove_extra_newlines=true`, `short_to_long_function_def=true`):
```bash
julia -e "using JuliaFormatter; format(\".\", verbose=true)"
```

### Documentation
```bash
julia --project=docs docs/make.jl
```

## Architecture

### Expression Tree Model

Users build expression trees from leaves (Variables, Constants) and nodes (Atoms). Each expression carries DCP metadata:

- **Vexity**: `ConstVexity`, `AffineVexity`, `ConvexVexity`, `ConcaveVexity`, `NotDcp` — composition rules in `src/dcp.jl`
- **Sign**: `Positive`, `Negative`, `NoSign`, `ComplexSign`
- **Monotonicity**: `Nondecreasing`, `Nonincreasing`, `ConstMonotonicity`, `NoMonotonicity`

A `Problem{T}` (`:minimize`, `:maximize`, `:satisfy`) wraps an objective expression and constraint list. DCP compliance is validated before solving.

### Solving Pipeline

When `solve!(problem, optimizer)` is called:

1. **Context creation** — `Context{T}` wraps an MOI model and tracks variable/constraint index mappings
2. **Recursive `conic_form!`** — traverses the expression tree, producing `SparseTape{T}` (lazy affine Ax+b representation) or `ComplexTape{T}` for complex expressions
3. **Reformulation** — atoms either decompose into simpler Convex primitives or directly add MOI conic constraints (SOC, SDP, exponential cone, etc.)
4. **MOI conversion** — `SparseTape` objects become `MOI.VectorAffineFunction` via `to_vaf()`
5. **Solution recovery** — solver results populate variable values and constraint duals through stored mappings

### Key Source Files

| File | Role |
|------|------|
| `src/dcp.jl` | DCP rule system (vexity/sign/monotonicity algebra) |
| `src/expressions.jl` | `AbstractExpr` base type |
| `src/variable.jl` | Variable types (`Variable`, `Semidefinite`, `ComplexVariable`, etc.) |
| `src/Constraint.jl` | `Constraint{S}` parameterized by MOI set type |
| `src/problems.jl` | `Problem{T}`, objective/constraint management, vexity checking |
| `src/SparseTape.jl` | Lazy affine transformation representation |
| `src/ComplexTape.jl` | Pair of SparseTapes for real/imaginary parts |
| `src/operate.jl`, `src/real_operate.jl`, `src/complex_operate.jl` | Tape arithmetic operations |
| `src/Context.jl` | `Context{T,M}` — solver context holding MOI model + metadata |
| `src/MOI_wrapper.jl` | `Optimizer{T,M}` — MOI-compatible wrapper |
| `src/solution.jl` | Solution extraction from solver |
| `src/atoms/` | 44 atom types (one file per atom) |
| `src/reformulations/` | Extended formulations (conv, dot, norm, quadform, etc.) |
| `src/sets/` | Custom MOI cone definitions |
| `src/problem_depot/` | Categorized test problems by cone type (affine, lp, socp, sdp, exp, mip) |

### Adding a New Atom

Each atom file in `src/atoms/` must implement:
- `sign(atom)` — output sign
- `monotonicity(atom)` — tuple of Monotonicity per child
- `curvature(atom)` — vexity of the atom itself
- `evaluate(atom)` — numeric evaluation on constants
- `new_conic_form!(context, atom)` — reformulation to conic constraints (the core method)

The `new_conic_form!` method either recursively builds from Convex primitives or directly adds MOI constraints to `context.model`. Use existing atoms like `NuclearNormAtom.jl` or `EuclideanNormAtom.jl` as templates.

Tests go in `test/test_atoms.jl` using the `_test_atom` and `_test_reformulation` helpers. Problem depot tests go in `src/problem_depot/problems/`.

### Compatibility

- Julia 1.6+
- Primary dependency: MathOptInterface v1.17
- CI runs on Julia 1.6, latest stable, and nightly
