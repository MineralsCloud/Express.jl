function color_string(str::AbstractString, color::AbstractChar)
    env = @match color begin
        'b' => "\033[34m"
        'r' => "\033[31m"
        'g' => "\033[32m"
        'm' => "\033[35m"
        _ => throw(ArgumentError("unknown color flag: $color"))
    end
    return env * str * "\033[0m\033[0m"
end

macro c_str(str, color = 'g')
    return color_string(str, first(color))
end # macro c_str
