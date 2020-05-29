module Workspaces

export LocalWorkspace, DockerWorkspace, SSHWorkspace

abstract type Workspace end
struct LocalWorkspace <: Workspace
    bin
    dir
    env
end
abstract type RemoteWorkspace <: Workspace end
abstract type VirtualMachineWorkspace <: RemoteWorkspace end
struct DockerWorkspace <: VirtualMachineWorkspace
    n
    container
end
struct SSHWorkspace <: RemoteWorkspace end

end # module Workspaces
