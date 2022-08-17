# MakeInput
jobify(f::UpdateTemplate, args...) = Job.(thunkify(f, args...))
