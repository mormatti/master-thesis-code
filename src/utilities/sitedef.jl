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