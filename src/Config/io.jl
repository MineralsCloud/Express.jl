using Configurations: @option
using Formatting: sprintf1

export DirStructure
export iofiles

@option struct InputFile
    name::String = "%s.in"
end

@option struct OutputFile
    name::String = "%s.out"
end

@option struct DirStructure
    rootdir::String = pwd()
    subdirname::String = "%s"
    input::InputFile = InputFile()
    output::OutputFile = OutputFile()
    DirStructure(rootdir, subdirname, input, output) =
        new(abspath(expanduser(rootdir)), subdirname, input, output)
end

function iofiles(ds::DirStructure, dirnames, filenames)
    dirs = map(dirnames) do dirname
        joinpath(ds.rootdir, sprintf1(ds.subdirname, dirname))
    end
    if length(filenames) == 1
        filename = fill(only(filenames), length(dirnames))
    end
    return map(dirs, filenames) do dir, filename
        input, output =
            sprintf1(ds.input.name, filename), sprintf1(ds.output.name, filename)
        joinpath(dir, input) => joinpath(dir, output)
    end
end
