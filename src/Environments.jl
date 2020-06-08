module Environments

export LocalEnvironment, DockerEnvironment, ServerEnvironment

abstract type CalculationEnvironment end
struct LocalEnvironment <: CalculationEnvironment
    n::Int
    bin::Any
    env::Any
end
abstract type RemoteEnvironment <: CalculationEnvironment end
abstract type VirtualMachineEnvironment <: RemoteEnvironment end
struct DockerEnvironment <: VirtualMachineEnvironment
    n::Int
    container::Any
    bin::Any
end
struct ServerEnvironment <: RemoteEnvironment end

end # module Environments
