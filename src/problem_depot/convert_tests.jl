function convert_test(file, prefix)
    # handle start and end:
    # delete @testset for... and final "end"
    file_contents = read(file, String)
    file_contents = replace(file_contents, r"\h*@testset.*for\ssolver.*\n" => "")
    file_contents = lstrip(file_contents)
    file_contents = file_contents[1:prevind(file_contents, first(findlast("end", file_contents)))]
    file_contents = replace(file_contents, r"\n\h\h\h\h(.*)" => SubstitutionString("\n\\1"))

    # handle special cases
    file_contents = special_cases(file_contents)

    # main work
    file_contents = replace_function_names(file_contents, prefix)
    file_contents = process_test_blocks(file_contents)
    file_contents = process_problems(file_contents)
    file_contents = replace(file_contents, "TOL" => "atol")

    open(joinpath("problems", "$(prefix).jl"), "w") do io
        write(io, file_contents)
    end
end

function get_indent(line)
    space_match = match(r"(\h*)\S", line)
    space_match === nothing ? "" : space_match.captures[]
end

function special_cases(file_contents)
    apply_each_line(file_contents) do line
        if strip(line) == "@test solve!(p, solver) === nothing"
            space = get_indent(line)
            return space * "output = solve!(p, solver)\n" * space * "@test output === nothing"
        elseif strip(line) == "@test_throws Exception solve!(p, solver)"
            space = get_indent(line)
            return space * "@test_throws Exception handle_problem!(p)"
        else
            return line
        end
    end
end

function replace_matched(f, str::AbstractString, regex::Regex, inner_regex::Regex = regex)
    g = str->begin
        m = match(inner_regex, str)
        if m === nothing
            @error "Match somehow not found" regex str
        end
        f(m)
    end
    replace(str, regex => g)
end

function replace_function_names(str, prefix)
    replace_matched(str, r"""@testset "(.*)" begin""") do m
        name = m.captures[]
        name = replace(name, " " => "_")
        name = replace(name, "#" => "")
        name = replace(name, "`" => "")
        output = "@add_problem $prefix function $(prefix)_$(name)(handle_problem!, valtest::Val{test} = Val(false), atol=1e-3, rtol=0.0, typ::Type{T} = Float64) where {T, test}"
    end
end

apply_each_line(f, str) = mapreduce(f, (x, y)->x * "\n" * y, split(str, '\n'); init = "")
function process_problems(str)
    apply_each_line(process_problem_line, str)
end

function process_test_blocks(str)
    indent = "    "
    replace_matched(str, r"(\n\h*?@test[\h(_throws)].*)+") do m
        space_match = match(r"\n(\h*)\S", m.match)
        if space_match === nothing
            space = ""
        else
            space = space_match.captures[]
        end
        test_statements = replace(m.match, r"\n(\h*\S)" => SubstitutionString("\n" * indent * "\\1"))
        test_statements = space * indent * lstrip(strip(test_statements))
        "\n" * space * "if test\n"  * test_statements * "\n" * space * "end"
    end
end

function process_problem_line(line)
    if occursin("@test", line)
        if !startswith(strip(line), "@test")
            @error "Line does not start with @test, but has @test" line
        end
        if occursin("solve!", line)
            @warn "Line has both `solve!` and `@test`; will need special handling" line
        end
        # new_line =  replace(line, "@test" => "test && @test")
        new_line = line
        if occursin("≈", new_line)
            if occursin(r"atol\h?=\h?TOL", new_line)
                new_line = replace(new_line, r"atol\h?=\h?TOL" => "atol=atol rtol=rtol")
            else
                new_line = new_line * " atol=atol rtol=rtol"
            end
        end
        return new_line
    elseif occursin("solve!", line)
        return replace_matched(line, r"solve!\((.+?)\,[^\)\,\n]*\)[^\)\,\n]*$") do m
            "handle_problem!($(m.captures[]))"
        end
    else
        return line
    end
end

for f in readdir("tests")
    if startswith(f, "test_")
        m = match(r"test_(.*)\.jl", f)
        prefix = m.captures[1]
        if prefix == "const"
            prefix = "constant"
        end
        convert_test(joinpath("tests", f), string(prefix))
    end
end
