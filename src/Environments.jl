module Environments

export LocalEnvironment, DockerEnvironment, ServerEnvironment

abstract type SimulationEnvironment end
struct LocalEnvironment <: SimulationEnvironment
    n::Int
    bin::Any
    dir::Any
    env::Any
end
abstract type RemoteEnvironment <: SimulationEnvironment end
abstract type VirtualMachineEnvironment <: RemoteEnvironment end
struct DockerEnvironment <: VirtualMachineEnvironment
    n::Int
    container::Any
end
struct ServerEnvironment <: RemoteEnvironment end

end # module Environments
