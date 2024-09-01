using ITensors

L = 100
g4 = 0.1
a = 1
g = sqrt(sqrt(g4))
sites = siteinds("Plaquette", L)

osH = opSumHamiltonian(N = L, g4 = g4)
H = MPO(osH, sites)

println("Computing DMRG...")
ψin = randomMPS(sites, L)
E0, ψ0 = dmrg(H, ψin;
            nsweeps = 10,
            maxdim = 100,
            cutoff = 10^(-13),
            outputlevel = 1,
)

localHamiltonians = []
println("Computing local MPOs...")
print("Done i = ")
for i in 1:L
    print("$i,")
    push!(localHamiltonians, MPO(opSumLocalHamiltonian(N = L, j = i, g4 = g4), sites))
end