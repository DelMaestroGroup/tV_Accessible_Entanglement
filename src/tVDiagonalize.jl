#Creates BoseHubbardDiagonalize module

module tVDiagonalize

using JeszenszkiBasis

export
    BdryCond,
    OBC,
    PBC,
    twositesRDM,
    TSRDM,
    NOTSRDM,
    NNNbipart,
    NNN,
    NONNN,

    sparse_hamiltonian,
    spatial_entropy,
    spatial_entropy_rhoij,
    spatial_entropyNNN

"""
Boundary conditions.
"""

@enum BdryCond PBC OBC
@doc "Periodic boundary conditions." PBC
@doc "Open boundary conditions." OBC

@enum twositesRDM NOTSRDM TSRDM
@doc "do not generate the two sites RDM." NOTSRDM
@doc "generate the two sites RDM." TSRDM

@enum NNNbipart NONNN NNN
@doc "do not include next nearest neighbor bipartition." NONNN
@doc "include next nearest neighbor bipartition." NNN


include("sparse_hamiltonian.jl")
include("spatial_entropy.jl")
include("spatial_entropy_rhoij.jl")
include("spatial_entropyNNN.jl")

end
