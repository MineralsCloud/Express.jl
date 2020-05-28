module Schemes

export LocalScheme, DockerScheme, SSHScheme

abstract type Scheme end
struct LocalScheme <: Scheme end
abstract type RemoteScheme <: Scheme end
abstract type VirtualMachineScheme <: RemoteScheme end
struct DockerScheme <: VirtualMachineScheme
    container::String
end
struct SSHScheme <: RemoteScheme end

end # module Schemes
