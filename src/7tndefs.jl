using ITensors

# We create a 2x2 complex matrix of zeros
id = zeros(Complex{Float64}, 2, 2)
id[1,1] = 1
id[2,2] = 1

sigmaPlus = zeros(Complex{Float64}, 2, 2)
sigmaPlus[1,2] = 1

Ul = 2 * sigmaPlus

sigmaZ = zeros(Complex{Float64}, 2, 2)
sigmaZ[1,1] = 1
sigmaZ[2,2] = -1

sigmaX = zeros(Complex{Float64}, 2, 2)
sigmaX[1,2] = 1
sigmaX[2,1] = 1

Etop = sigmaZ / 2

# We define the site types corresponding to the spin 3/2
ITensors.space(::SiteType"Plaquette") = 2
ITensors.op(::OpName"I",::SiteType"Plaquette") = id
ITensors.op(::OpName"U",::SiteType"Plaquette") = Ul
ITensors.op(::OpName"U",::SiteType"Plaquette") = Ul
ITensors.op(::OpName"sigma+",::SiteType"Plaquette") = sigmaPlus
ITensors.op(::OpName"Udag",::SiteType"Plaquette") = Ul'
ITensors.op(::OpName"sigma-",::SiteType"Plaquette") = sigmaPlus'
ITensors.op(::OpName"Etop",::SiteType"Plaquette") = Etop
ITensors.op(::OpName"sigmaZ",::SiteType"Plaquette") = sigmaZ
ITensors.op(::OpName"U+Udag",::SiteType"Plaquette") = Ul + Ul'
ITensors.op(::OpName"sigmaX",::SiteType"Plaquette") = sigmaX
ITensors.op(::OpName"Wl",::SiteType"Plaquette") = sigmaPlus - sigmaPlus'

function opSumHamiltonianDirichlet(; L::Integer, a::Real = 1, g4::Real, TwoETop0::Integer, TwoETopLp1::Integer)
    g2 = sqrt(g4)
    os = OpSum()

    # Electric term
    os += (g2/(2*a)),"I",1
    os += -TwoETop0 * (g2/(2*a)),"Etop",1
    for j=1:L-1
        os += (g2/(2*a)),"I",j
        os += -2 * (g2/(2*a)),"Etop",j,"Etop",j+1
    end
    os += (g2/(2*a)),"I",L
    os += -TwoETopLp1 * (g2/(2*a)),"Etop",L
    os += -1/2 * (g2/(2*a)),"I",L

    # magnetic term
    for j=1:L
        os += -(1/(2*g2*a)),"U+Udag",j
    end
    
    return os
end

function opSumLocalHamiltonianDirichlet(; L::Integer, a::Real = 1, g4::Real, TwoETop0::Integer, TwoETopLp1::Integer, j::Integer)
    g2 = sqrt(g4)
    os = OpSum()
    if j == 1
        os += 5/4*(g2/(2*a)),"I",1
        os += -(g2/(2*a)),"Etop",1,"Etop",2
        os += -TwoETop0*(g2/(2*a)),"Etop",1
        os += -(1/(2*g2*a)),"U+Udag",1
    elseif j == L
        os += 5/4*(g2/(2*a)),"I",L
        os += -(g2/(2*a)),"Etop",L,"Etop",L-1
        os += -TwoETopLp1*(g2/(2*a)),"Etop",L
        os += -(1/(2*g2*a)),"U+Udag",L
    else
        os += (g2/(2*a)),"I",j
        os += -(g2/(2*a)),"Etop",j-1,"Etop",j
        os += -(g2/(2*a)),"Etop",j,"Etop",j+1
        os += -(1/(2*g2*a)),"U+Udag",j
    end
    return os
end

function cartesianProduct(vectOfVects1, vectOfVects2)
    result = []
    for i in eachindex(vectOfVects1)
        for j in eachindex(vectOfVects2)
            push!(result, vcat(vectOfVects2[j], vectOfVects1[i]))
        end
    end
    return result
end

function cartesianPower(vectOfVects, n)
    result = deepcopy(vectOfVects)
    for i in 2:n
        result = cartesianProduct(deepcopy(result), deepcopy(vectOfVects))
    end
    return result
end

function rawWannierLocalized(; j::Integer)
    os = OpSum()
    os += 1,"Wl",j
    return os
end

function wannierLocalized(; ext::Integer, j::Real, coefficientsArray)
    ext = abs(ext)
    os = OpSum()
    opStringsArray = cartesianPower([["I"], ["sigma+"]], 2 * ext + 1)
    @assert length(coefficientsArray) == length(opStringsArray)
    for jw in eachindex(opStringsArray)
        opt = opStringsArray[jw]
        if ext == 0
            os += coefficientsArray[jw],opt[1],j
        elseif ext == 1
            os += coefficientsArray[jw],opt[1],j-1,opt[2],j,opt[3],j+1
        elseif ext == 2
            os += coefficientsArray[jw],opt[1],j-2,opt[2],j-1,opt[3],j,opt[4],j+1,opt[5],j+2
        elseif ext == 3
            os += coefficientsArray[jw],opt[1],j-3,opt[2],j-2,opt[3],j-1,opt[4],j,opt[5],j+1,opt[6],j+2,opt[7],j+3
        elseif ext == 4
            os += coefficientsArray[jw],opt[1],j-4,opt[2],j-3,opt[3],j-2,opt[4],j-1,opt[5],j,opt[6],j+1,opt[7],j+2,opt[8],j+3,opt[9],j+4
        elseif ext == 5
            os += coefficientsArray[jw],opt[1],j-5,opt[2],j-4,opt[3],j-3,opt[4],j-2,opt[5],j-1,opt[6],j,opt[7],j+1,opt[8],j+2,opt[9],j+3,opt[10],j+4,opt[11],j+5
        elseif ext == 6
            os += coefficientsArray[jw],opt[1],j-6,opt[2],j-5,opt[3],j-4,opt[4],j-3,opt[5],j-2,opt[6],j-1,opt[7],j,opt[8],j+1,opt[9],j+2,opt[10],j+3,opt[11],j+4,opt[12],j+5,opt[13],j+6
        else
            error("Not implemented.")
        end
    end
    return os
end

function rawGaussianPacket(; L::Integer, jbar::Integer, pbar::Real, sigma::Real)
    a(j) = exp(im * pbar * j) * exp(-(j - jbar)^2 / (2 * sigma^2))
    os = OpSum()
    for j=1:L
        os += a(j),"Wl",j
    end
    return os
end

function gaussianPacket(; L::Integer, ext::Integer, jbar::Real, pbar::Real, sigma::Real, coefficientsArray)
    ext = abs(ext)
    a(j) = exp(im * pbar * j) * exp(-(j - jbar)^2 / (2 * sigma^2))
    os = OpSum()
    opStringsArray = cartesianPower([["I"], ["sigma+"]], 2 * ext + 1)
    @assert length(coefficientsArray) == length(opStringsArray)
    for j=1+ext:L-ext
        for jw in eachindex(opStringsArray)
            opt = opStringsArray[jw]
            if ext == 0
                os += a(j)*coefficientsArray[jw],opt[1],j
            elseif ext == 1
                os += a(j)*coefficientsArray[jw],opt[1],j-1,opt[2],j,opt[3],j+1
            elseif ext == 2
                os += a(j)*coefficientsArray[jw],opt[1],j-2,opt[2],j-1,opt[3],j,opt[4],j+1,opt[5],j+2
            elseif ext == 3
                os += a(j)*coefficientsArray[jw],opt[1],j-3,opt[2],j-2,opt[3],j-1,opt[4],j,opt[5],j+1,opt[6],j+2,opt[7],j+3
            elseif ext == 4
                os += a(j)*coefficientsArray[jw],opt[1],j-4,opt[2],j-3,opt[3],j-2,opt[4],j-1,opt[5],j,opt[6],j+1,opt[7],j+2,opt[8],j+3,opt[9],j+4
            elseif ext == 5
                os += a(j)*coefficientsArray[jw],opt[1],j-5,opt[2],j-4,opt[3],j-3,opt[4],j-2,opt[5],j-1,opt[6],j,opt[7],j+1,opt[8],j+2,opt[9],j+3,opt[10],j+4,opt[11],j+5
            elseif ext == 6
                os += a(j)*coefficientsArray[jw],opt[1],j-6,opt[2],j-5,opt[3],j-4,opt[4],j-3,opt[5],j-2,opt[6],j-1,opt[7],j,opt[8],j+1,opt[9],j+2,opt[10],j+3,opt[11],j+4,opt[12],j+5,opt[13],j+6
            else
                error("Not implemented.")
            end
        end
    end 
    return os
end


function photonicPacket(; L::Integer, ext::Integer, jbar::Real, pmin::Real, pmax::Real, coefficientsArray)
    ext = abs(ext)
    sign = Int(round(pmax / abs(pmax)))
    Nmin = Int(ceil(L/2 * (abs(pmin) / π)))
    Nmax = Int(ceil(L/2 * (abs(pmax) / π)))
    a(j) = sum([exp(im*2*π*n*sign/L * (j - jbar)) for n in Nmin:Nmax])
    os = OpSum()
    opStringsArray = cartesianPower([["I"], ["sigma+"]], 2 * ext + 1)
    @assert length(coefficientsArray) == length(opStringsArray)
    for j=1+ext:L-ext
        for jw in eachindex(opStringsArray)
            opt = opStringsArray[jw]
            if ext == 0
                os += a(j)*coefficientsArray[jw],opt[1],j
            elseif ext == 1
                os += a(j)*coefficientsArray[jw],opt[1],j-1,opt[2],j,opt[3],j+1
            elseif ext == 2
                os += a(j)*coefficientsArray[jw],opt[1],j-2,opt[2],j-1,opt[3],j,opt[4],j+1,opt[5],j+2
            elseif ext == 3
                os += a(j)*coefficientsArray[jw],opt[1],j-3,opt[2],j-2,opt[3],j-1,opt[4],j,opt[5],j+1,opt[6],j+2,opt[7],j+3
            elseif ext == 4
                os += a(j)*coefficientsArray[jw],opt[1],j-4,opt[2],j-3,opt[3],j-2,opt[4],j-1,opt[5],j,opt[6],j+1,opt[7],j+2,opt[8],j+3,opt[9],j+4
            elseif ext == 5
                os += a(j)*coefficientsArray[jw],opt[1],j-5,opt[2],j-4,opt[3],j-3,opt[4],j-2,opt[5],j-1,opt[6],j,opt[7],j+1,opt[8],j+2,opt[9],j+3,opt[10],j+4,opt[11],j+5
            elseif ext == 6
                os += a(j)*coefficientsArray[jw],opt[1],j-6,opt[2],j-5,opt[3],j-4,opt[4],j-3,opt[5],j-2,opt[6],j-1,opt[7],j,opt[8],j+1,opt[9],j+2,opt[10],j+3,opt[11],j+4,opt[12],j+5,opt[13],j+6
            else
                error("Not implemented.")
            end
        end
    end 
    return os
end