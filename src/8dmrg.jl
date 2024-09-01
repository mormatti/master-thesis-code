using ITensors
using JLD2

a = 1

LArrayTot = [100, 300]
g4ArrayTot = [0.0001, 0.001, 0.01, 0.1, 0.2]
BCArrayTot = ["dirichlet1", "dirichlet2", "vonNeumann1", "vonNeumann2"]

LArray = [100, 300]
g4Array = [0.0001, 0.001, 0.01, 0.1, 0.2]
BCArray = ["dirichlet1"]

supersuperpath = "/Users/Mattia/Code/Quantools.jl/scattering/photon-photon/data/pbc/nosector/dsh1dsv2"

for L in LArray

    localPath = "$supersuperpath/L$L"

    # Create folder of L if it doesn't exist
    if !isdir(localPath)
        println("L = $L folder not found. Creating it...")
        mkdir(localPath)
    else
        println("L = $L folder found.")
    end

    # Create sites file if it doesn't exist
    if !isfile("$localPath/sites.jld2")
        println("Sites file not found. Creating it...")
        sites = siteinds("Plaquette", L)
        save_object("$localPath/sites.jld2", sites)
    else
        println("Sites file found.")
        sites = load_object("$localPath/sites.jld2")
    end

    jHalf = Int(floor(L/2))

    for g4 in g4Array

        localPath = "$supersuperpath/L$L/g4=$g4"

        # if the folder of g4 doesn't exist, create it
        if !isdir(localPath)
            println("g4 = $g4 folder not found. Creating it...")
            mkdir(localPath)
        else
            println("g4 = $g4 folder found.")
        end

        for BC in BCArray

            localPath = "$supersuperpath/L$L/g4=$g4/$BC"

            # if the folder of BC doesn't exist, create it
            if !isdir(localPath)
                println("$BC folder not found. Creating it...")
                mkdir(localPath)
            else
                println("$BC folder found.")
            end

            # Create MPO file if it doesn't exist
            if !isfile("$localPath/hamiltonianMPO.jld2")
                println("MPO file not found. Creating it...")
                if BC == "dirichlet1"
                    OpSumH = opSumHamiltonianDirichlet(L = L, g4 = g4, TwoETop0 = 1, TwoETopLp1 = 1)
                elseif BC == "dirichlet2"
                    OpSumH = opSumHamiltonianDirichlet(L = L, g4 = g4, TwoETop0 = 1, TwoETopLp1 = -1)
                end
                hamiltonianMPO = MPO(OpSumH, sites)
                save_object("$localPath/hamiltonianMPO.jld2", hamiltonianMPO)
            else
                println("MPO file found.")
                hamiltonianMPO = load_object("$localPath/hamiltonianMPO.jld2")
            end

            # Searching the ground state with DMRG if it doesn't exist
            if !isfile("$localPath/groundState.jld2")
                println("Ground state file not found. Searching for it...")
                ψin = randomMPS(sites, L)
                E0, ψ0 = dmrg(hamiltonianMPO, ψin;
                    nsweeps = 50,
                    maxdim = 200,
                    cutoff = 10^(-13),
                    outputlevel = 1,
                )
                save_object("$localPath/groundState.jld2", ψ0)
            else
                println("Ground state file found.")
                ψ0 = load_object("$localPath/groundState.jld2")
            end
        end # BC loop
    end # g4 loop
end # L loop