using LinearAlgebra

spin12Set = [-1/2, 1/2]
spin1Set = [-1, 0, 1]
spin32Set = [-3/2, -1/2, 1/2, 3/2]
spin2Set = [-2, -1, 0, 1, 2]
spin52Set = [-5/2, -3/2, -1/2, 1/2, 3/2, 5/2]
spin3Set = [-3, -2, -1, 0, 1, 2, 3]
matterSet = [0,1]
antimatterSet = [-1,0]


sVSet = spin1Set

sT1Set = sVSet
sB1Set = sVSet
ψT1Set = matterSet
ψB1Set = antimatterSet

sT2Set = sVSet
sB2Set = sVSet
ψT2Set = antimatterSet
ψB2Set = matterSet

dim1 = length(sT1Set) * length(sB1Set) * length(ψT1Set) * length(ψB1Set)
dim2 = length(sT2Set) * length(sB2Set) * length(ψT2Set) * length(ψB2Set)

MA = zeros(Int, dim1, dim2)
MB = zeros(Int, dim1, dim2)

i1 = i2 = 0

for sT1 in sT1Set
    for sB1 in sB1Set
        for ψT1 in ψT1Set
            for ψB1 in ψB1Set
                global i1 = i1 + 1
                for sT2 in sT2Set
                    for sB2 in sB2Set
                        for ψT2 in ψT2Set
                            for ψB2 in ψB2Set
                                global i2 = i2 + 1
                                if sT2 + sB2 - (sT1 + sB1) == ψT2 + ψB2
                                    MA[i1,i2] = 1
                                end
                            end
                        end
                    end
                end
                println("")
                global i2 = 0
            end
        end
    end
end
i2 = 0



i1 = i2 = 0

for sT1 in sT2Set
    for sB1 in sB2Set
        for ψT1 in ψT2Set
            for ψB1 in ψB2Set
                global i1 = i1 + 1
                for sT2 in sT1Set
                    for sB2 in sB1Set
                        for ψT2 in ψT1Set
                            for ψB2 in ψB1Set
                                global i2 = i2 + 1
                                if sT2 + sB2 - (sT1 + sB1) == ψT2 + ψB2
                                    MB[i1,i2] = 1
                                end
                            end
                        end
                    end
                end
                println("")
                global i2 = 0
            end
        end
    end
end
i2 = 0

for val in eigen(MA).values
    if abs(val) > 10^(-15)
        println(val)
    end
end