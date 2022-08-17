using Configurations: @option
using Formatting: sprintf1

@option struct InputFile
    name::String = "%s.in"
end

@option struct OutputFile
    name::String = "%s.out"
end

@option struct SubDir
    root::String = pwd()
    name::String = "%s"
    input::InputFile = InputFile()
    output::OutputFile = OutputFile()
end

function getfiles(dir::SubDir, name, filename)
    path = joinpath(dir.root, sprintf1(dir.name, name))
    input, output = sprintf1(dir.input.name, filename), sprintf1(dir.output.name, filename)
    return joinpath(path, input) => joinpath(path, output)
end
