using ITensors
using Plots
using LaTeXStrings
using JLD2

ACTIVATE_PROCESS = true

if ACTIVATE_PROCESS

LArray = [100] # 100 standard
g4Array = [0.1] # 0.1 standard
sigmaArray = [3] # 3 standard
momentumOverπArray = [1/2] # 1/2 standard
scattering = true
a = 1

for L in LArray
for g4 in g4Array
for sigma in sigmaArray
for momentumOverPi in momentumOverπArray

# Parameters
g = sqrt(sqrt(g4))
# Boundary conditions
BC = "dirichlet1"
twoEtop0 = 1
twoEtopLp1 = 1
# Wannier parameters
ext = 2
w = 2 * ext + 1
# Gaussian packet parameters
jbar = 20
momentumAbs = momentumOverPi * π
# pmin = π/4
# pmax = π/2
# TEBD parameters
cutoff = 1E-10
dt = 0.05 # 0.05 standard
ttotal = 400 # 400
maxdim = 50
timeSteps = 0:dt:ttotal
timeResolution = 300 # 300 standard
# Plotting variables
timeGroupingNumber = Int(ceil(length(timeSteps) / timeResolution))
numberOfGroups = Int(ceil(length(timeSteps) / timeGroupingNumber))
plotEnergies = []
plotEnergiesLog = []

# Loading objects
print("Loading objects...")
supersuperpath = "/Users/Mattia/Code/Quantools.jl/scattering/photon-photon/data/pbc/nosector/dsh1dsv2"
sites = load_object("$supersuperpath/L$L/sites.jld2")
ψ0 = load_object("$supersuperpath/L$L/g4=$g4/$BC/groundState.jld2")
HMPO = load_object("$supersuperpath/L$L/g4=$g4/$BC/hamiltonianMPO.jld2")
coeffArray = load_object("$supersuperpath/L13/g4=$g4/wannierCoefficients.jld2")
println("Done.")

# Computing MPOs
print("Creating MPOs...")
MPO1 = MPO(gaussianPacket(; L = L, ext = ext, jbar = jbar, pbar = momentumAbs, sigma = sigma, coefficientsArray = coeffArray), sites)
MPO2 = MPO(gaussianPacket(; L = L, ext = ext, jbar = L-jbar, pbar = -momentumAbs, sigma = sigma, coefficientsArray = coeffArray), sites)
# MPO1 = MPO(rawGaussianPacket(; L = L, jbar = jbar, pbar = momentumAbs, sigma = sigma), sites)
# MPO2 = MPO(rawGaussianPacket(; L = L, jbar = L-jbar, pbar = -momentumAbs, sigma = sigma), sites)
# MPO1 = MPO(wannierLocalized(; ext = ext, j = jbar, coefficientsArray = coeffArray), sites, splitblocks=false)
# MPO1 = MPO(rawWannierLocalized(; j = Integer(L/2)), sites)
# MPO1 = MPO(photonicPacket(; L = L, ext = ext, jbar = jbar, pmin = pmin, pmax = pmax, coefficientsArray = coeffArray), sites)
# display(MPO1)
println("Done.")

# Applying MPOs to groundstate
print("Applying MPOs to groundstate...")
ψ1 = apply(MPO1, ψ0, maxdim = maxdim)
ψ1 = apply(MPO2, ψ1, maxdim = maxdim)
ψ1 = ψ1 / norm(ψ1)
println("Done.")

# Computing local Hamiltonians
print("Computing local Hamiltonians...")
localHamiltonians = []
println("Computing local MPOs...")
print("Step: 1/$L ")
for i in 1:L
    print("\r")
    print("Step: $i/$L ")
    push!(localHamiltonians, MPO(opSumLocalHamiltonianDirichlet(; L = L, g4 = g4, TwoETop0 = 1, TwoETopLp1 = 1, j = i), sites))
end
println("Done.")

# Plotting energy densities
print("Plotting the energy densities of the initial state...")
plot(xtickfontsize = 14,
ytickfontsize = 14,
dpi = 300)
E1minusE0 = real(inner(ψ1, HMPO, ψ1) - inner(ψ0, HMPO, ψ0))
energyDensities = []
println("Computing energy densies...")
print("Step: 1/$L ")
for i=1:L
    print("\r")
    print("Step: $i/$L ")
    push!(energyDensities, abs((real(inner(ψ1, localHamiltonians[i], ψ1) - inner(ψ0, localHamiltonians[i], ψ0)))/(E1minusE0)))
end
plot!(1:L,
    energyDensities,
    marker=(:circle,3))
savefig("EnergyDensityPlot.png")

plot!(yaxis=:log)

savefig("EnergyDensityLog.png")

println("Done.")

ψ2 = deepcopy(ψ1)

# Computing all the TEBD gates
print("Computing all the TEBD gates...")
gates = ITensor[]
I1 = op("I", sites[1])
I2 = op("I", sites[2])
E1 = op("Etop", sites[1])
E2 = op("Etop", sites[2])
U1 = op("U+Udag", sites[1])
U2 = op("U+Udag", sites[2])
hj = g^2/(2*a) * 2 * (7/8 * I1 * I2 - twoEtop0/2 * E1 * I2 - E1 * E2) - 1/(2*a*g^2) * U1 * I2
Gj = exp(-im * dt / 2 * hj)
push!(gates, Gj)
for j in 2:(L-2)
    local s1 = sites[j]
    local s2 = sites[j+1]
    local Ij = op("I", s1)
    local Ij1 = op("I", s2)
    local Ej = op("Etop", s1)
    local Ej1 = op("Etop", s2)
    local Uj = op("U+Udag", s1)
    local Uj1 = op("U+Udag", s2)
    local hj = g^2/(2*a) * (Ij * Ij1 - 2 * Ej * Ej1) - 1/(2*a*g^2) * Uj * Ij1
    local Gj = exp(-im * dt / 2 * hj)
    push!(gates, Gj)
end
ILm1 = op("I", sites[L-1])
IL = op("I", sites[L])
ELm1 = op("Etop", sites[L-1])
EL = op("Etop", sites[L])
ULm1 = op("U+Udag", sites[L-1])
UL = op("U+Udag", sites[L])
hj = g^2/(2*a) * 2 * (7/8 * ILm1 * IL - ELm1 * EL - ILm1 * EL * twoEtopLp1/2) - 1/(2*a*g^2) * ULm1 * IL - 1/(2*a*g^2) * ILm1 * UL
Gj = exp(-im * dt / 2 * hj)
push!(gates, Gj)
append!(gates, reverse(gates)) # Include gates in reverse order too
println("Done.")


function convert_from_seconds(sec::Int)
    x, seconds = divrem(sec, 60)
    hours, minutes = divrem(x, 60)
    hours, minutes, seconds
end

println("Computing evolution...")
# Compute and print <Sz> at each time step
# then apply the gates to go to the next time

print("number of groups = $numberOfGroups")

timeIndex = 1
print("Estimated remaining time: ")
for groupIndex in 1:numberOfGroups
    Δt = @elapsed begin
        for groupTimeIndex in 1:timeGroupingNumber
            timeIndex > length(timeSteps) && break
            t = timeSteps[timeIndex]
            ψ2 = apply(gates, ψ2; cutoff, maxdim)
            normalize!(ψ2)
            timeIndex += 1
        end
        energyDensities = []
        energyDensitiesLog = []
        for sitepos in 1:L
            energyDensity = real(inner(ψ2, localHamiltonians[sitepos], ψ2) - inner(ψ0, localHamiltonians[sitepos], ψ0))
            push!(energyDensities, energyDensity)
            push!(energyDensitiesLog, log10(abs(energyDensity)))
        end
        push!(plotEnergies, energyDensities)
        push!(plotEnergiesLog, energyDensitiesLog)
    end
    ttotal = Δt * (numberOfGroups - groupIndex)
    ts = convert_from_seconds(Int(round(ttotal)))
    print("\r")
    print("Estimated remaining time: $(ts[1]):$(ts[2]):$(ts[3]).")
    # We save the plot energies array
    supersuperpath = "/Users/Mattia/Code/Quantools.jl/scattering/photon-photon/data/pbc/nosector/dsh1dsv2"
    save_object("$supersuperpath/L$L/g4=$g4/dirichlet1/plotEnergies,σ=$sigma,j=$jbar,kovπ=$momentumOverPi,scat=$scattering.jld2", plotEnergies)
end

end # momentum 
end # sigma
end # g4
end # L
end # ACTIVATE_PROCESS