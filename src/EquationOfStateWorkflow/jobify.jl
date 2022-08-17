# MakeInput
jobify(f::MakeInput, args...) = Job.(thunkify(f, args...))
