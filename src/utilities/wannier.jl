using Optim
using LinearAlgebra
using LaTeXStrings


function momentumEnergy(z::Complex{Float64}, L::Int)
    r, phi = abs(z), angle(z)
    n = phi / (pi / L)
    n = round(Int, n)
    if n % 2 == 1 || n % 2 == -1
        if n >= 0
            n = n - L
        else
            n = n + L
        end
        r = -r
    end
    n = n / 2
    return r, n
end


function momentum(z::Complex{Float64}, L::Int)
    return momentumEnergy(z, L)[2]
end


function energy(z::Complex{Float64}, L::Int)
    return momentumEnergy(z, L)[1]
end



# Let n be a positive integer. The followin function, given a (2n+1)x(2n+1) matrix and
# an array of n real numbers, it returns...
function functional(Λ::Matrix, θ::Vector)
    return real(sum(Λ[i, j] * exp(im * (θ[j] - θ[i])) for i in 1:size(θ, 1), j in 1:size(θ, 1)))
end

# The following function doubles the array and symmetrizes it by adding the elements
# in reverse order. For example, [1, 2] becomes [2, 1, 0, 1, 2].
function doubling_symmetrize(a::Vector)
    return [reverse(a);[0];a]
end


# The following function returns the value of the functional for the symmetric case.
function symmetricFunctional(Λ::Matrix, θ::Vector)
    return functional(Λ, doubling_symmetrize(θ))
end


function optimumThetas(Λ::Matrix)
    f(θ) = symmetricFunctional(Λ, θ)
    n = Int(floor(size(Λ, 1)/2))
    θ = rand(Float64, n)
    result = optimize(f, θ, NelderMead())
    return result.minimizer
end


function evolve(ψ, t, L, ħ = 1)
    F = load_object("data/pbc/nosector/dsh1dsv2/L$L/eigenHT.jld2")
    energies = load_object("data/pbc/nosector/dsh1dsv2/L$L/energies.jld2")
    evolutionMatrix = F.vectors
    diagonalMatrix = diagm(energies)
    return evolutionMatrix * exp(-im * (1 / ħ) * t * diagonalMatrix) * evolutionMatrix' * ψ
end


# The following function takes as input a list of vectors. Each vector represents a plot.
# Each vector is a list of points. For example, the vector [1,2,3] is a plot with three points.
# The function returns the gif with the plot of the points, each vector of the list represents
# a frame of the gif animation.
function gifPlot(listOfVectors::Vector, filename::String, yMin::Real, yMax::Real)
    anim = @animate for i ∈ 1:length(listOfVectors)
        plot(listOfVectors[i], ylim = (yMin,yMax), legend = false, xlabel = "plaquette site", ylabel = "energy density")
    end
    gif(anim, "anim_fps15.gif", fps = 20)
end




#=

ACTIVATE_PROCESS = true

if ACTIVATE_PROCESS

plot_font = "Computer Modern"
default(fontfamily=plot_font, linewidth=2, framestyle=:box, label=nothing, grid=false)

recomputeHamiltonian = false
recomputeLocalHamiltonian = true
recumputeHT = false
recomputeEigenHT = false

dsh = 1
dsv = 2

LArray = [9]
gArray = [0.0001, 0.001, 0.01, 0.1, 0.2, 0.5, 1.0]

for L in LArray
for g4 in gArray
@time begin

local g = sqrt(sqrt(g4))
local CB = load_object("data/pbc/nosector/dsh$(dsh)dsv$(dsv)/L$(L)/chainBasis.jld2")
local dirPath = "data/pbc/nosector/dsh$(dsh)dsv$(dsv)/L$(L)/g4=$(g4)"

# Directory
if !isdir(dirPath)
    print("Creating folder... ")
    mkdir(dirPath)
    println("Done.")
end

# Hamiltonian operator
local H = []
if !isfile("$(dirPath)/hamiltonian.jld2") || recomputeHamiltonian
    print("Constructing Hamiltonian operator... ")
    H = total_hamiltonian(CB, g)
    save_object("$(dirPath)/hamiltonian.jld2", H)
else
    print("Loading Hamiltonian operator... ")
    H = load_object("$(dirPath)/hamiltonian.jld2")
end
println("Done.")

# Central local Hamiltonian
if !isfile("$(dirPath)/centralLocalhamiltonian.jld2") || recomputeLocalHamiltonian
    print("Constructing the central local Hamiltonian... ")
    Hc = local_total_hamiltonian(CB, g, Int((L+1)/2))
    save_object("$(dirPath)/centralLocalhamiltonian.jld2", Hc)
else
    print("Loading the central local Hamiltonian... ")
    Hc = load_object("$(dirPath)/centralLocalhamiltonian.jld2")
end
println("Done.")

# Symmetry operators
print("Loading symmetry operators... ")
local T = load_object("data/pbc/nosector/dsh$(dsh)dsv$(dsv)/L$(L)/translationOperator.jld2")
local R = load_object("data/pbc/nosector/dsh$(dsh)dsv$(dsv)/L$(L)/reflectionOperator.jld2")
println("Done.")

# Hamiltonian times translation operator
if !isfile("$(dirPath)/HT.jld2") || recumputeHT
    print("Constructing Hamiltonian times translation operator... ")
    HT = H * T
    save_object("$(dirPath)/HT.jld2", HT)
else
    print("Loading Hamiltonian times translation operator... ")
    HT = load_object("$(dirPath)/HT.jld2")
end
println("Done.")

# Diagonalizing Hamiltonian times translation operator
if !isfile("$(dirPath)/eigenHT.jld2") || recomputeEigenHT
    print("Diagonalizing Hamiltonian times translation operator")
    eigenHT = eigen(HT)
    save_object("$(dirPath)/eigenHT.jld2", eigenHT)
else
    print("Loading eigensystem of H * T...")
    eigenHT = load_object("$(dirPath)/eigenHT.jld2")
end
eigenHTvalues = eigenHT.values
eigenHTvectors = eigenHT.vectors
println("Done.")

# Computing energies and momenta
print("Computing energies and momenta... ")
energies = []
momenta = []
for i in 1:length(eigenHTvalues)
    push!(energies, energy(eigenHTvalues[i], L))
    push!(momenta, momentum(eigenHTvalues[i], L))
end
save_object("$(dirPath)/energies.jld2", energies)
save_object("$(dirPath)/momenta.jld2", momenta)
kmax = (L-1)/2
kmin = -kmax
println("Done.")

# Computing groundstate energy and state
print("Computing groundstate energy and state... ")
groundStateIndex = findmin(energies)[2]
maximumEnergyIndex = findmax(energies)[2]
groundStateEnergy = energies[groundStateIndex]
save_object("$(dirPath)/groundStateEnergy.jld2", groundStateEnergy)
groundStateVector = eigenHTvectors[:, groundStateIndex]
save_object("$(dirPath)/groundStateVector.jld2", groundStateVector)
maximumEnergy = energies[maximumEnergyIndex]
save_object("$(dirPath)/maximumEnergy.jld2", maximumEnergy)
println("Done.")

# Computing first band momenta, energies and states
print("Computing first band momenta, energies and states... ")
firstBandMomenta = kmin:kmax
firstBandIndices = []
firstBandEnergies = []
firstBandStates = []
for k in kmin:kmax
    indicesK = findall(x -> x == k, momenta)
    arr = energies[indicesK]
    minIndex = findmin(arr)[2]
    if k == 0
        energiesK = energies[indicesK]
        minIndex = findfirst(x -> x == minimum(filter(y -> y != minimum(arr), arr)), arr)
    end
    firstBandIndices = [firstBandIndices; indicesK[minIndex]]
    push!(firstBandEnergies, energies[indicesK[minIndex]])
    push!(firstBandStates, eigenHTvectors[:,indicesK[minIndex]])
end
save_object("$(dirPath)/firstBandMomenta.jld2", firstBandMomenta)
save_object("$(dirPath)/firstBandEnergies.jld2", firstBandEnergies)
save_object("$(dirPath)/firstBandStates.jld2", firstBandStates)
println("Done.")

print("Plotting first band... ")
plot(firstBandMomenta, firstBandEnergies, legend = false, xlabel = "momentum", ylabel = "energy", title = "First band energies")
savefig("data/pbc/nosector/dsh$(dsh)dsv$(dsv)/L$(L)/g4=$(g4)/firstBandEnergies.png")
println("Done.")

# Computing the symmetrized first band states
print("Computing the symmetrized first band states... ")
symmetrizedFirstBandStates = []
for i in 1:L
    ψ = firstBandStates[i]
    ψm = firstBandStates[L + 1 - i]
    phase = (i <= Int((L+1)/2)) ? 1.0 : -1.0 / (ψm' * R * ψ)
    symmetrizedψ = phase * ψ
    push!(symmetrizedFirstBandStates, symmetrizedψ)
end
save_object("$(dirPath)/symmetrizedFirstBandStates.jld2", symmetrizedFirstBandStates)
println("Done.")

# Computing the Λ matrix
print("Computing the Λ matrix... ")
Λ = zeros(Complex, length(firstBandMomenta), length(firstBandMomenta))
Γ(x) = 1.0/L * sum([j^2 * exp(2 * im * j * π * x / L) for j in kmin:kmax])
for i1 in eachindex(firstBandMomenta)
    for i2 in eachindex(firstBandMomenta)
        k1 = firstBandMomenta[i1]
        k2 = firstBandMomenta[i2]
        Λ[i1,i2] = (symmetrizedFirstBandStates[i1]' * Hc * symmetrizedFirstBandStates[i2]) * Γ(k1 - k2)
    end
end
println("Done.")

# Computing the maximally localized Wannier function
print("Computing the Wannier... ")
θ = doubling_symmetrize(optimumThetas(Λ))
local Wjc = sum([((1/√L) * exp(im * θ[i]) * symmetrizedFirstBandStates[i]) for i in eachindex(θ)])
save_object("$(dirPath)/centralWannier.jld2", Wjc)
println("Done.")

# Computing the energy density of the Wannier
print("Computing the energy density of the Wannier... ")
energydensity = [real(Wjc' * Hc * Wjc) - groundStateEnergy/L]
global ψA = deepcopy(Wjc)
global ψB = deepcopy(Wjc)
for j in 1:((L-1)/2)
    ψA = T * ψA
    ψB = T' * ψB
    energydensity = [real(ψB' * Hc * ψB) - groundStateEnergy/L; energydensity; real(ψA' * Hc * ψA) - groundStateEnergy/L]
end
save_object("$(dirPath)/energydensityWannier.jld2", energydensity)
energydensityNormalized = abs.(energydensity / (sum(energydensity)))
save_object("$(dirPath)/energydensityNormalizedWannier.jld2", energydensityNormalized)
plot(energydensityNormalized,
    yaxis = :log,
    legend = false,
    xlabel = L"plaquette site $j$",
    ylabel = (L"\langle\hat H_j \rangle_W - \langle \hat H_j \rangle_\Omega / (\langle \hat H \rangle_W - \langle \hat H \rangle_\Omega)"), 
    title = "Energy density of the Wannier", dpi = 300)
savefig("$(dirPath)/energydensity.png")
println("Done.")

end #time
end #g for cycle
end #L for cycle


# Plotting the energy density of the Wannier
LArray = [13]
plot(xlabel = L"Plaquette site $j$",
    yaxis = :log,
    ylabel = L"(\langle \hat H_j \rangle_W - \langle \hat H_j \rangle_\Omega) / ( \langle \hat H \rangle_W - \langle \hat H \rangle_\Omega)", 
    title = "Normalized energy density of the Wannier",
    dpi = 300,
    xticks = 1:15)
for L in LArray
    for g4 in gArray
        energyDensity = load_object("data/pbc/nosector/dsh1dsv2/L$(L)/g4=$(g4)/energydensityNormalizedWannier.jld2")
        plot!(energyDensity, label = latexstring("g^4 = $(g4)"))
    end
end
savefig("wannier.png")

end # ACTIVATE_PROCESS

=#
