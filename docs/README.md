# Documentation and examples

This file describes how to generate Convex.jl's documentation and examples. If
you just want to read the documentation, or try out the examples, please visit
the website: <https://www.juliaopt.org/Convex.jl/stable/>. You can download a
zip file of Jupyter notebooks of all the examples there as well.

The examples and documentation are constructed together. Run

```sh
julia --project=docs -e 'using Pkg; Pkg.instantiate(); include("docs/make.jl")'
```

to generate the examples notebooks (which will be placed in `docs/notebooks`)
and the documentation itself, which is generated into the `doc/build` folder,
and can be previewed by opening a webserver there. Note that this command can
take some time. To generate the documentation without updating the examples,
set `ENV["CONVEX_SKIP_EXAMPLES"]="true"` before including `docs/make.jl`.

To generate a single Jupyter notebook, run e.g.

```julia
Literate.notebook(file_path, notebook_dir, execute=false) # or execute = true, to run the code
```

Then the notebook can be opened with IJulia.

To just generate a single markdown file, run e.g.

```julia
fix_math_md(content) = replace(content, r"\$\$(.*?)\$\$"s => s"```math\1```")
Literate.markdown(file_path, output_directory; preprocess = fix_math_md)
```

This won't execute the code however; that is done by Documenter, so the above
`docs/make.jl` file is needed. This can be slow because all the examples will be
re-run. By re-including `docs/make.jl` into a running session, however, compile
times can be minimized.

The `fix_math_md` function allows us to use `$$` for LaTeX display in the
Literate.jl files by replacing it with the tags expected by Documenter.jl when
we generate the markdown.
