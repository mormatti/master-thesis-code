using DelimitedFiles, DataFrames
using Plots
using LaTeXStrings
using JLD2

plot_font = "Computer Modern"
default(fontfamily=plot_font, linewidth=1, label=nothing, background_color=RGBA(1,1,1,0))

cutoff = 1E-10
dt = 0.05 # 0.05 standard
ttotal = 400 # 400
maxdim = 50
ttotal = 400 # 400
timeSteps = 0:dt:ttotal
timeResolution = 300 # 300 standard
# Plotting variables
timeGroupingNumber = Int(ceil(length(timeSteps) / timeResolution))
numberOfGroups = Int(ceil(length(timeSteps) / timeGroupingNumber))

timeArray = [0, ttotal/4, ttotal/2, 3*ttotal/4, ttotal]
timeArrayN = timeArray / ttotal * numberOfGroups

# system params
L = 100 # 100 standard
g4Array = [0.1] # 0.1 standard
sigmaArray = [3] # 3 standard
kovπArray = [0.5] # 1/2 standard
jbar = 20 # 20 standard
scattering = true

for g4 in g4Array
for sigma in sigmaArray
for kovπ in kovπArray

supersuperpath = "/Users/Mattia/Code/Quantools.jl/scattering/photon-photon/data/pbc/nosector/dsh1dsv2"
plotEnergies = load_object("$supersuperpath/L$L/g4=$g4/dirichlet1/plotEnergies,σ=$sigma,j=$jbar,kovπ=$kovπ,scat=$scattering.jld2")

plotEnergies = plotEnergies / maximum(plotEnergies[1])
fps = 30 # Int(ceil(timeResolution / 15))
display(plotEnergies)
oneRealSecondIs = 1 / (ttotal / numberOfGroups / timeGroupingNumber)
anim = @animate for animIndex in eachindex(plotEnergies)
    plot(plotEnergies[animIndex],
        ylim = (-0.1,2.2), legend = false,
        dpi = 600, grid = false, ratio = 5,
        framestyle = :none, axis = false
        )
    print("$animIndex ")
end
print(" Finished!")
gif(anim, "animationg4=$g4,sigma=$sigma,kovπ=$kovπ,scattering=$scattering.gif", fps = fps)

#= # Now we transform the plotEnergies array into a matrix
# and plot it
plotEnergies = hcat(plotEnergies...)
heatmap(plotEnergies, 
        yflip = true, 
        xticks=(timeArrayN, timeArray),
        xtickfontsize = 14,
        ytickfontsize = 14,
        dpi = 300)
savefig("EnergyDensityPlotHeatmap.png")

# We consider the log of the energy densities
plotEnergiesLog = hcat(plotEnergiesLog...)
heatmap(plotEnergiesLog,
        yflip = true,
        xticks=(timeArrayN, timeArray),
        xtickfontsize = 14,
        ytickfontsize = 14,
        dpi = 300)
savefig("EnergyDensityPlotHeatmapLog.png") =#

end # kovπ
end # sigma
end # g4


#= # we create an array of strings like "data/exactDiagonalization/dsh1/nosector/LN/data.txt"
# where N ranging from 2 to 8

# filepaths = ["data/DMRG/dsh1/L$N/data16.txt" for N in [10]]
filepaths = ["data/DMRG/dsh4/L10/bd300.txt", "data/DMRG/dsh4/L50/bd300.txt"]

# We create an array of 10 colors, which forms a transition from black to blue
# colors = [RGB(0.3 + 0.7 * 1/(length(filepaths) - i + 1), 0, 0) for i in 1:length(filepaths)]
colors = [RGB(0, 0, 0), RGB(1, 0, 0)]
# colors = [RGB(0.8,0,0.2), RGB(0,0.8,0.2), RGB(0,0.2,0.8)]

plot(xlabel = L"\lambda", ylabel = L"\partial_\lambda^2 E_0 / (g^2/4a) / L", label = "Energy")

# We set the y range from -100 to 100
# plot!(xlim = (0, 1))
# plot!(ylim = (-30, 10))

# We erase the legend from the plot
plot!(legend = false)

plot!(size=(800,500))
plot!(xguidefontsize=15, yguidefontsize=15, tickfontsize=8)

for i in eachindex(filepaths)
    local myarray = open(readdlm,filepaths[i])
    local λarray = myarray[:,1]
    local Earray = myarray[:,2]

    local n = 2

    λarray = λarray[1:n:end]
    Earray = Earray[1:n:end]

    local sep = λarray[2] - λarray[1]

    # We compute the second derivative of the energy array
    Earray = discreteForwardSecondDerivative(Earray, sep)
    λarray = eliminateFirstAndLastElement(λarray)

    #We plot! λarray, energyarray changing the color of the line from blue to lightblue
    plot!(λarray, Earray, color = colors[i])
end

savefig("Energy.pdf") =#