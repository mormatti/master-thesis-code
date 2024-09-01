using Optim
using LinearAlgebra
using LinearSolve

# The following function, given a set S of vectors, check how many are linearly independent
function linearlyIndependentElements(S)
    # we form the matrix of the vectors in S
    M = zeros(Complex{Float64}, length(S[1]), length(S))
    for i in 1:length(S)
        M[:, i] = S[i]
    end
    return rank(M)
end

# The following function shows the eigenvalues of the matrix formed by the scalar products
# of the vectors in S
function scalarProductMatrix(S)
    spm = zeros(Complex{Float64}, length(S), length(S))
    for i in 1:length(S)
        for j in 1:length(S)
            spm[i, j] = S[i]' * S[j]
        end
    end
    return spm
end

function scalarProductVector(S, v)
    spv = zeros(Complex{Float64}, length(S))
    for i in 1:length(S)
        spv[i] = v' * S[i]
    end
    return spv
end

function interpolate(vector, VSet)
    M = scalarProductMatrix(VSet)
    v = scalarProductVector(VSet, vector)
    coefficients = inv(M) * v
    vInterpolated = zeros(Complex{Float64}, length(vector))
    for i in 1:length(VSet)
        vInterpolated += coefficients[i] * VSet[i]
    end
    println("Norm of the difference: $(norm(vector - vInterpolated))")
    println("scalar product: $(vector' * vInterpolated)")
    return coefficients, vInterpolated
end

function vectorDistance1(v1, v2)
    return norm(v1 - v2)
end

function vectorDistance2(v1, v2)
    return norm(sqrt(1 - (abs(v1' * v2))^2 / (norm(v1)^2 * norm(v2)^2)))
end

function applyUnitary(n::Int, U, O, v)
    w = deepcopy(v)
    if n > 0
        for i in 1:n
            w = U * w
        end
        w = O * w
        for i in 1:n
            w = U' * w
        end
    elseif n < 0
        for i in 1:n
            w = U' * w
        end
        w = O * w
        for i in 1:n
            w = U * w
        end
    end
    return w
end

include("5wannier.jl")

#=
dsh, dsv = 1, 2
L = 11 # An odd number
ext = 3 # The extension of the interpolating operator
w = 2 * ext + 1 # An odd number
g4 = 0.0001
jc = Integer(floor(L / 2) + 1)
=#

dsh, dsv = 1, 2

LArray =  [13]
for Li in eachindex(LArray)
L = LArray[Li]
println("L = $L")
jc = Integer(floor(L / 2) + 1)

supersuperpath = "/Users/Mattia/Code/Quantools.jl/scattering/photon-photon/data/pbc/nosector/dsh$(dsh)dsv$(dsv)"
dirPathUp = "$supersuperpath/L$L"
CB = load_object("$(dirPathUp)/chainBasis.jld2")
T = load_object("$(dirPathUp)/translationOperator.jld2")
R = load_object("$(dirPathUp)/reflectionOperator.jld2")

a(j) = 1/2 * click_operator(CB, j)
b(j) = 1/2 * clock_operator(CB, j)
ab(j) = a(j) * b(j)
A = []
for i in 1:L
    push!(A, a(i))
end

g4Array = [0.0001, 0.001, 0.01, 0.1, 0.2]
for g4i in eachindex(g4Array)
g4 = g4Array[g4i] 
println("g4 = $g4")


dirPath = "$supersuperpath/L$L/g4=$g4"
Wjc = load_object("$(dirPath)/centralWannier.jld2")
Ω = load_object("$(dirPath)/groundStateVector.jld2")
F = load_object("$(dirPath)/eigenHT.jld2")
Emin = load_object("$(dirPath)/groundStateEnergy.jld2")
Emax = load_object("$(dirPath)/maximumEnergy.jld2")
H = load_object("$(dirPath)/hamiltonian.jld2")
V = F.vectors
en = load_object("$(dirPath)/energies.jld2")
val = F.values
dimH = length(Ω)

extArray = [2]
for ext in extArray
w = 2 * ext + 1 # An odd number
println("w = $w")

vectorSet = [Ω]
for i in (-ext):(ext)
    toAdd = []
    for v in vectorSet
        # w = applyUnitary(i, T, a(jc), v)
        push!(toAdd, A[jc+i] * v)
    end
    for v in toAdd
        push!(vectorSet, v)
    end
end
# println("Number of vectors in vectorSet: $(length(vectorSet))")

# println("Computing coefficients...")
# println(eigen(scalarProductMatrix(vectorSet)).values)
coeff, iterp = interpolate(Wjc, vectorSet)

# println(coeff)

# Results[Li, g4i, exti, 1] = vectorDistance1(Wjc, iterp)
# Results[Li, g4i, exti, 2] = vectorDistance2(Wjc, iterp)

# print("Vector dist.1: $(vectorDistance1((-I + 2 * a(jc)) * Ω, Wjc)), ")
# println("Vector dist.2: $(vectorDistance2((-I + 2 * a(jc)) * Ω, Wjc)) ")

# We save coeff in a file
save_object("$(dirPath)/wannierCoefficients.jld2", coeff)

end # for ext
end # for L
end # for g4

# save_object("data/pbc/nosector/dsh$(dsh)dsv$(dsv)/OperatorResults.jld2", Results)

# display(Results)