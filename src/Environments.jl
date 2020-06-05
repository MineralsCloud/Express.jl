module Environments

export LocalEnvironment, DockerEnvironment, ServerEnvironment

abstract type SimulationEnvironment end
struct LocalEnvironment <: SimulationEnvironment
    n::Int
    bin::Any
    env::Any
end
abstract type RemoteEnvironment <: SimulationEnvironment end
abstract type VirtualMachineEnvironment <: RemoteEnvironment end
struct DockerEnvironment <: VirtualMachineEnvironment
    n::Int
    container::Any
    bin::Any
end
struct ServerEnvironment <: RemoteEnvironment end

end # module Environments
