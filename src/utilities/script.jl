using LaTeXStrings
using Plots
using JLD2

# Write here your script code.

# We save in a string the path of the file OperatorResultsSaved.jld2
dsh = 1
dsv = 2
path = "data/pbc/nosector/dsh$(dsh)dsv$(dsv)/OperatorResultsSaved.jld2"

# We load the file OperatorResultsSaved.jld2
ResultsLoaded = load_object(path)

LArray =  [3, 5, 7, 9, 11, 13]
g4Array = [0.0001, 0.001, 0.01, 0.1, 0.2, 0.5, 1.0]
extArray = [0, 1, 2, 3, 4, 5, 6]

plot_font = "Computer Modern"
default(fontfamily=plot_font, linewidth=2, framestyle=:box, label=nothing, grid=false)

# We set the y ticks (in log) to be 10^(-1), 10^(-2), ..., 10^(-5)
plot(xlabel = L"Extension $w$", ylabel = L"Interpolation error $\sqrt{1 - \langle W | W' \rangle / |W|^2 |W'|^2 }$", yaxis = :log10, dpi = 300)
plot!(xticks = (1:13, 1:13))
plot!(yticks = 10.0 .^ collect(range(-11,-1,step=1)))

for Li in [6]
    L = LArray[Li]
    for g4i in 1:7
        g4 = g4Array[g4i]
        # We square all the elements of ResultsLoaded[Li, g4i, :, 2]
        ResultsLoadedSquared = ResultsLoaded[Li, g4i, :, 2] .^ 2
        plot!(1:2:13, ResultsLoadedSquared, label = latexstring("g^4 = $g4"), title = latexstring("L = $L"))
    end
end
savefig("data/pbc/nosector/dsh$(dsh)dsv$(dsv)/OperatorResultsPlot.png")
