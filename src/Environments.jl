module Workspaces

export LocalWorkspace, DockerWorkspace, SSHWorkspace

abstract type Workspace end
struct LocalWorkspace <: Workspace
    n::Int
    bin
    dir
    env
end
abstract type RemoteWorkspace <: Workspace end
abstract type VirtualMachineWorkspace <: RemoteWorkspace end
struct DockerWorkspace <: VirtualMachineWorkspace
    n::Int
    container
end
struct SSHWorkspace <: RemoteWorkspace end

end # module Workspaces
