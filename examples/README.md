# Examples

These examples are available from the documentation website:
<https://www.juliaopt.org/Convex.jl/stable/examples/>.

The examples are written as
[Literate.jl](https://github.com/fredrikekre/Literate.jl) files. They can be
built into Jupyter notebooks by running

```julia
julia --project=. build.jl
```

from this directory. This will (re-)generate `.ipynb` files in the `notebooks`
directory. To add a new example or update one of the
examples, please edit `literate` files and submit a pull request. The notebooks
will be built automatically when the documentation is generated and added to
<https://www.juliaopt.org/Convex.jl/dev/examples/>.
