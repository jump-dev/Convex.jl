# Documentation and examples

The examples and documentation are constructed together. Run

```sh
julia --project=docs -e 'using Pkg; Pkg.instantiate(); include("docs/make.jl")'
```

to generate the examples notebooks (which will be placed in `docs/notebooks`) and the documentation itself, which is generated into the `doc/build` folder, and can be previewed by opening a webserver there.
