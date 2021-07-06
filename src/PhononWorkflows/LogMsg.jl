@action LogMsg
function (x::LogMsg{T})(start = true) where {T}
    startend = start ? "starts" : "ends"
    with_logger(current_logger()) do
        println("The calculation $T $startend at $(format(now(), "HH:MM:SS u dd, yyyy")).")
    end
end

buildjob(x::LogMsg, start) = InternalAtomicJob(() -> x(start))
