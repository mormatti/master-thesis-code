"""This file contains the definition of the chainBasis object and the functions that
operate on it. A chainBasis is the set of all possible chains that can be formed by
concatenating plaquettes with the same spin representation."""

include("0spinZ.jl")
include("1plaquette.jl")
include("2chain.jl")


using JLD2

mutable struct ChainBasis
    chains::Vector{Chain}

    """The constructor for chainBasis. It takes a vector of chains."""
    ChainBasis(chains::Vector{Chain}) = new(chains)
end

"""It returns a copy of the chainBasis"""
function Base.deepcopy(cb::ChainBasis)::ChainBasis
    chains = Chain[]
    for c in cb.chains
        push!(chains, deepcopy(c))
    end
    return ChainBasis(chains)
end

"""It returns the number of chains in the chainBasis"""
Base.length(cb::ChainBasis) = length(cb.chains)

"""The constructor for chainBasis. It takes no arguments and returns the void chainBasis."""
ChainBasis() = ChainBasis(Chain[])

"""given a spin representation with two integer for the horizontal and vertical doubled spin,
it returns the chainBasis formed of all the possible chain with single plaquettes with that
spin representation"""
function plaquetteBasis(dsh::Int, dsv::Int)::ChainBasis
    chains = Chain[]
    for p in Plaquette_list(dsh, dsv, dsh, dsv)
        push!(chains, Chain(p))
    end
    return ChainBasis(chains)
end

"""given a spin representation with two integer for the horizontal and vertical doubled spin,
it returns the chainBasis formed of all the possible chain with semi-plaquettes with that
spin representation"""
function semiPlaquetteBasis(dsh::Int, dsv::Int)::ChainBasis
    chains = Chain[]
    for p in semiPlaquette_list(dsh, dsv, dsh, dsv)
        push!(chains, Chain(p))
    end
    return ChainBasis(chains)
end

"""Prints in a readable way all the chains in the chainBasis"""
function Base.show(cb::ChainBasis; io::IO = stdout)
    println(io, "chainBasis:")
    for c in cb.chains
        println(io, " ")
        show(c, io = io)
    end
end

"""It returns the chainBasis formed by all the possible gauge concatenations of the chains in
the two chainBasis. If two chains are not compatible, no new chain is added to the chainBasis"""
function ⊗(cb1::ChainBasis, cb2::ChainBasis)::ChainBasis
    chains = Chain[]
    for c1 in cb1.chains
        for c2 in cb2.chains
            if c1 ▷ c2
                if c1 * c2 in chains
                    continue
                else
                    push!(chains, c1 * c2)
                end 
            end
        end
    end
    return ChainBasis(chains)
end

function clever_gauge_concatenation(cb1::ChainBasis, cb2::ChainBasis)::ChainBasis
    chains = Chain[]
    for c1 in cb1.chains
        for c2 in cb2.chains
            if c1 ≥ c2
                push!(chains, c1 ⊕ c2)
            end
        end
    end
    return ChainBasis(chains)
end

"""It eliminates from the chainBasis all the chains that are not under periodic boundary
conditions"""
function gauge_closure(cb::ChainBasis)::ChainBasis
    cbNew = ChainBasis(filter(x -> gauge_lockable(x), cb.chains))
    for c in cbNew.chains
        lockChain(c)
    end
    return cbNew
end

function clever_gauge_closure(cb::ChainBasis)::ChainBasis
    cbNew = deepcopy(cb)

    # We create a void array to store indices
    indices = Int[]

    for i in eachindex(cbNew.chains)
        if gauge_lockable(cbNew.chains[i])
            gauge_lock(cbNew.chains[i])
        else
            # we add the index i to the array indices
            push!(indices, i)
        end
    end

    # we delete all the chains that are not gauge lockable
    deleteat!(cbNew.chains, indices)

    return cbNew
end

"""This function eliminates all the chains that do not belongs to a certain given polarization
sector, i.e. all the plaquettes has to have the same given longitudinal polarization"""
function select_polarization_sector(cb::ChainBasis, sector::Int)::ChainBasis
    return ChainBasis(filter(chain -> belongs_to_longpol_sector(chain, sector), cb.chains))
end

"""It returns a chiainBasis with the periodic boundary conditions. It takes as arguments 
the horizontal and vertical doubled spin, the length of the chain the polarization sector 
and the selection of the polarization sector"""
function compute_chainBasis(; 
    horizontalDoubleSpin::Int = 1,
    verticalDoubleSpin::Int = 2,
    chainLength::Int = 5,
    polarizationSector = nothing
    )::ChainBasis

    PB = semiPlaquetteBasis(horizontalDoubleSpin, verticalDoubleSpin)
    CB = deepcopy(PB)
    for i in 1:chainLength - 1
        CB = clever_gauge_concatenation(CB,PB)
        if polarizationSector !== nothing
            CB = select_polarization_sector(CB, polarizationSector)
        end
    end
    CB = clever_gauge_closure(CB)
    return CB
end

#= dsh = 1
dsv = 2
L = 13
dirPathUp = "/Users/Mattia/Code/Quantools.jl/scattering/photon-photon/data/pbc/nosector/dsh$(dsh)dsv$(dsv)/L$L"
CB = compute_chainBasis(horizontalDoubleSpin = 1, verticalDoubleSpin = 2, chainLength = 13)
println(length(CB))
save_object("$dirPathUp/chainBasis.jld2", CB) =#