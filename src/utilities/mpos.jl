using ITensors

function opSumHamiltonian(; N::Int64, g4 = 1.0, a = 1)
    g2 = sqrt(g4)
    os = OpSum()
    for j=2:N-1
        os += (g2/(2*a)),"I",j
        os += -(g2/(2*a)),"Etop",j-1,"Etop",j
        os += -(g2/(2*a)),"Etop",j,"Etop",j+1
    end
    for j=1:N
        os += -(1/(2*g2*a)),"U+Udag",j
    end
    return os
end

function opSumLocalHamiltonian(; N::Int64, j::Int64 = 1, g4 = 1.0, a = 1)
    g2 = sqrt(g4)
    os = OpSum()
    if j == 1
        os += -(g2/(2*a)),"Etop",1,"Etop",2
        os += -(1/(2*g2*a)),"U+Udag",1
    elseif j == N
        os += -(g2/(2*a)),"Etop",N-1,"Etop",N
        os += -(1/(2*g2*a)),"U+Udag",N
    else
        os += (g2/(2*a)),"I",j
        os += -(g2/(2*a)),"Etop",j-1,"Etop",j
        os += -(g2/(2*a)),"Etop",j,"Etop",j+1
        os += -(1/(2*g2*a)),"U+Udag",j
    end
    return os
end

function opSumWl(; N::Int64, j::Int64 = 1)
    os = OpSum()
    os += 1,"Wl",j
    return os
end

function rawPacket(; N::Int64, jbar::Int64 = 1, pbar::Real = 1.0, sigma::Real = 1.0)
    os = OpSum()
    a(j) = exp(im * pbar * j) * exp(-(j - jbar)^2 / (2 * sigma^2))
    for j=1:N
        os += a(j),"Wl",j
    end
    return os
end