using Literate

literate_path = joinpath(@__DIR__(), "literate")
if !(@isdefined(notebooks_path))
    notebooks_path = joinpath(@__DIR__(), "notebooks")
end

rm(notebooks_path, recursive=true, force=true)
mkdir(notebooks_path)

for dir in readdir(literate_path)
    dir_path = joinpath(literate_path, dir)
    isdir(dir_path) || continue
    @info "Processing directory $dir"
    notebook_dir = joinpath(notebooks_path, dir)
    isdir(notebook_dir) || mkdir(notebook_dir)
    for file in readdir(dir_path)
        file_path = joinpath(dir_path, file)
            out_path = joinpath(notebooks_path, dir, file)
        if endswith(file, ".jl")
            Literate.notebook(file_path, notebook_dir, execute=false)
        else
            cp(file_path, out_path)
        end
    end
end
