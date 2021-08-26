module Recipes

function buildworkflow(cfgfile)
    config = loadconfig(cfgfile)
    if isfile(config.recover)
        w = deserialize(config.recover)
        typeassert(w, Workflow)
        return w
    else
        settings = load(cfgfile)
        x = if settings["workflow"] == "phonon dispersion"
            PhononDispersion
        elseif settings["workflow"] == "vdos"
            VDos
        else
            error("unsupported option!")
        end
        return begin
            @intjob(LogMsg{Scf}()(true)) ▷ @intjob(MakeInput{Scf}()(cfgfile)) ▷
            buildjob(MakeCmd{Scf}(), cfgfile) ▷ @intjob(LogMsg{Scf}()(false)) ▷
            @intjob(LogMsg{Dfpt}()(true)) ▷ @intjob(MakeInput{Dfpt}()(cfgfile)) ▷
            buildjob(MakeCmd{Dfpt}(), cfgfile) ▷ @intjob(LogMsg{Dfpt}()(false)) ▷
            @intjob(LogMsg{RealSpaceForceConstants}()(true)) ▷
            @intjob(MakeInput{RealSpaceForceConstants}()(cfgfile)) ▷
            buildjob(MakeCmd{RealSpaceForceConstants}(), cfgfile) ▷
            @intjob(LogMsg{RealSpaceForceConstants}()(false)) ▷ @intjob(LogMsg{x}()(true)) ▷
            @intjob(MakeInput{x}()(cfgfile)) ▷ buildjob(MakeCmd{x}(), cfgfile) ▷
            @intjob(LogMsg{x}()(false))
        end
    end
end

end
