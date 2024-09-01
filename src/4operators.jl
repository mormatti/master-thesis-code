using LinearAlgebra

include("3chainBasis.jl")

using JLD2
using Plots

"""It returns the electric energy of a chain. It is the sum of the squares of the spin 
z projections of all the link of the chain. This method takes into account the periodicity
of the chain."""
function electric_energy(chain::Chain)::Float64
    if chain == Chain()
        return 0.0
    else
        energy = 0.0
        for p in chain.plaquettes
            energy += (p.sT.dsz / 2)^2 + (p.sL.dsz / 2)^2 + (p.sB.dsz / 2)^2
        end
        if !chain.closed
            energy += (last(chain).sR.dsz / 2)^2
        end
        return energy
    end
end

function local_electric_energy(chain::Chain, j::Integer)
    if chain == Chain()
        return 0.0
    else
        p = chain.plaquettes[j]
        return 0.5 * (p.sL.dsz / 2)^2 + (p.sT.dsz / 2)^2 + (p.sB.dsz / 2)^2 + 0.5 * (p.sR.dsz / 2)^2
    end
end

"""Given a chain basis, it returns the electric Hamiltonian term of the system."""
function electric_hamiltonian_term(basis::ChainBasis)::Matrix{Float64}
    H = zeros(Float64, length(basis), length(basis))
    for i in eachindex(basis.chains)
        H[i, i] = electric_energy(basis.chains[i])
    end
    return H
end

function local_electric_hamiltonian_term(basis::ChainBasis, j::Integer)::Matrix{Float64}
    H = zeros(Float64, length(basis), length(basis))
    for i in eachindex(basis.chains)
        H[i, i] = local_electric_energy(basis.chains[i], j)
    end
    return H
end


"""Given a chain basis, it returns the magnetic Hamiltonian term of the system."""
function magnetic_hamiltonian_term(basis::ChainBasis)::Matrix{Float64}
    H = zeros(Float64, length(basis), length(basis))
    C = basis.chains
    for i in eachindex(C)
        chain = C[i]
        for k in eachindex(chain.plaquettes)
            if U_applicable(chain, k)
                chainU = deepcopy(chain)
                U(chainU, k)
                j = findfirst(==(chainU), C)
                if j !== nothing
                    H[i, j] = U_prefactor(chain, k)
                end
            end
            if Udag_applicable(chain, k)
                chainUdag = deepcopy(chain)
                Udag(chainUdag, k)
                j = findfirst(==(chainUdag), C)
                if j !== nothing
                    H[i, j] = Udag_prefactor(chain, k)
                end
            end
        end
    end
    return H
end


function local_magnetic_hamiltonian_term(basis::ChainBasis, k::Integer)::Matrix{Float64}
    H = zeros(Float64, length(basis), length(basis))
    C = basis.chains
    for i in eachindex(C)
        chain = C[i]
        if U_applicable(chain, k)
            chainU = deepcopy(chain)
            U(chainU, k)
            j = findfirst(==(chainU), C)
            if j !== nothing
                H[i, j] = U_prefactor(chain, k)
            end
        end
        if Udag_applicable(chain, k)
            chainUdag = deepcopy(chain)
            Udag(chainUdag, k)
            j = findfirst(==(chainUdag), C)
            if j !== nothing
                H[i, j] = Udag_prefactor(chain, k)
            end
        end
    end
    return H
end

function click_operator(basis::ChainBasis, k::Integer)::Matrix{Float64}
    click = zeros(Float64, length(basis), length(basis))
    C = basis.chains
    for i in eachindex(C)
        chain = C[i]
        if U_applicable(chain, k)
            chainU = deepcopy(chain)
            U(chainU, k)
            j = findfirst(==(chainU), C)
            if j !== nothing
                click[i, j] = U_prefactor(chain, k)
            end
        end
    end
    return click
end

function clock_operator(basis::ChainBasis, k::Integer)::Matrix{Float64}
    clock = zeros(Float64, length(basis), length(basis))
    C = basis.chains
    for i in eachindex(C)
        chain = C[i]
        if Udag_applicable(chain, k)
            chainUdag = deepcopy(chain)
            Udag(chainUdag, k)
            j = findfirst(==(chainUdag), C)
            if j !== nothing
                clock[i, j] = Udag_prefactor(chain, k)
            end
        end
    end
    return clock
end


"""Given a chain basis, it returns the translation operator of the system."""
function translation_operator(basis::ChainBasis)::Matrix{Float64}
    H = zeros(Float64, length(basis), length(basis))
    C = basis.chains
    for i in eachindex(C)
        chain = C[i]
        if chain.closed
            chainT = translated(chain)
            j = findfirst(==(chainT), C)
            if j !== nothing
                H[i, j] = 1.0
            end
        end
    end
    return H
end


"""Given a chain basis, it returns the reflection operator of the system."""
function reflection_operator(basis::ChainBasis)::Matrix{Float64}
    H = zeros(Float64, length(basis), length(basis))
    C = basis.chains
    for i in eachindex(C)
        chain = C[i]
        if chain.closed
            chainR = y_reflected(chain)
            j = findfirst(==(chainR), C)
            if j !== nothing
                H[i, j] = 1.0
            end
        end
    end
    return H
end


function total_hamiltonian(basis::ChainBasis, g::Float64)::Matrix{Float64}
    println("Ok minus!")
    return g^2 / 2 * electric_hamiltonian_term(basis) - 1 / (2 * g^2) * magnetic_hamiltonian_term(basis)
end


function local_total_hamiltonian(basis::ChainBasis, g::Float64, j::Integer)::Matrix{Float64}
    println("Ok minus!")
    return g^2 / 2 * local_electric_hamiltonian_term(basis, j) - 1 / (2 * g^2) * local_magnetic_hamiltonian_term(basis, j)
end

#= L = 13
dsh, dsv = 1, 2
superpath = "/Users/Mattia/Code/Quantools.jl/scattering/photon-photon/data/pbc/nosector/dsh$(dsh)dsv$(dsv)/L$(L)"
CB = load_object("$superpath/chainBasis.jld2")
T = translation_operator(CB)
R = reflection_operator(CB)
save_object("$superpath/translationOperator.jld2", T)
save_object("$superpath/reflectionOperator.jld2", R) =#