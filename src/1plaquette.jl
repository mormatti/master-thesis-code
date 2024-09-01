"""The plaquette type. It has four fields, one for each side of the plaquette. The sides are
ordered B, R, T, L, where B is the bottom side, R is the right side, T is the top side, and
L is the left side."""

include("0spinZ.jl")

mutable struct Plaquette
    sB::SpinZ
    sR::SpinZ
    sT::SpinZ
    sL::SpinZ

    """The constructor for plaquette. It takes the spins in the order B, R, T, L."""
    function Plaquette(sB::SpinZ, sR::SpinZ, sT::SpinZ, sL::SpinZ)
        new(sB, sR, sT, sL)
    end
end

"""The longitudinal polarization is the sum of the top and bottom horizontal spin z projections
of the plaquette. It is a symmetry invariant quantity in a chain of plaquettes."""
longitudinal_polarization(plaquette::Plaquette) = plaquette.sB.dsz + plaquette.sT.dsz

"""Returns a deep-copy of the plaquette."""
function Base.deepcopy(p::Plaquette)::Plaquette
    return Plaquette(deepcopy(p.sB), deepcopy(p.sR), deepcopy(p.sT), deepcopy(p.sL))
end

"""The equality operator for plaquette. It checks if the spins are equal for both plaquettes."""
function Base.:(==)(p1::Plaquette, p2::Plaquette)::Bool
    p1.sB == p2.sB && p1.sR == p2.sR && p1.sT == p2.sT && p1.sL == p2.sL
end

"""Given a spin representation (4 values representing the spins of the four sides of the plaquette,
in the order B, R, T, L), it returns a random plaquette with the given spin representation,
where the random projection spins are chosen uniformly at random."""
function Plaquette_random(dsB::Int, dsR::Int, dsT::Int, dsL::Int)::Plaquette
    sBr = SpinZ_random(dsB)
    sRr = SpinZ_random(dsR)
    sTr = SpinZ_random(dsT)
    sLr = SpinZ_random(dsL)
    return Plaquette(sBr, sRr, sTr, sLr)
end

"""Returns a random plaquette given the spin representations of the horizontal and 
vertical spins of the plaquette."""
Plaquette_random(dsh::Int, dsv::Int)::Plaquette = Plaquette_random(dsh, dsv, dsh, dsv)

"""The "precedes" binary operator, it returns true if the common rung between the two plaquettes
has the same spin"""
Base.:>(p1::Plaquette, p2::Plaquette)::Bool = p1.sR.dsz == p2.sL.dsz

"""The "gauge-precedes" binary operator, it returns true if the common rung between the two plaquettes
has the same spin, and the nodes are consistent with the Gauss' law"""
▷(p1::Plaquette, p2::Plaquette)::Bool = p2.sT.dsz - p1.sT.dsz == p1.sB.dsz - p2.sB.dsz == p1.sR.dsz == p2.sL.dsz

"""The "weak precedes" binary operator. It returns true if the longitudinal polarization of the
two plaquettes are the same, and the common rung has a spin value which is consistent with its maximum spin."""
Base.:≥(p1::Plaquette, p2::Plaquette)::Bool = 
    (p1.sR.ds == p2.sL.ds
    && longitudinal_polarization(p1) == longitudinal_polarization(p2)
    && abs(p1.sT.dsz - p2.sT.dsz) <= p1.sR.ds)

"""Given the four values of the spins of the four sides of the plaquette, it returns
a list of all the different plaquettes with all the combinations of spin projections.
The order of the plaquettes in the list is the same as the order of the spins in the input."""
function Plaquette_list(dsB::Int, dsR::Int, dsT::Int, dsL::Int)::Vector{Plaquette}
    plaquettes = []
    for sB in SpinZ_list(dsB)
        for sR in SpinZ_list(dsR)
            for sT in SpinZ_list(dsT)
                for sL in SpinZ_list(dsL)
                    push!(plaquettes, Plaquette(sB, sR, sT, sL))
                end
            end
        end
    end
    return plaquettes
end


"""Given the four values of the spins of the four sides of the plaquette, it returns
a list of all the different plaquettes with all the combinations of the top and bottom
horizontal spin projections. The left and right spins projections are taken random."""
function semiPlaquette_list(dsB::Int, dsR::Int, dsT::Int, dsL::Int)::Vector{Plaquette}
    plaquettes = []
    for sB in SpinZ_list(dsB)
        for sT in SpinZ_list(dsT)
            sL = SpinZ_random(dsL)
            sR = SpinZ_random(dsR)
            push!(plaquettes, Plaquette(sB, sR, sT, sL))
        end
    end
    return plaquettes
end


"""Returns true if the plaquette operator U can be applied to the plaquette, i.e. if the 
plaquette spins are incrementable / decrementable in a loop."""
function U_applicable(p::Plaquette)::Bool
    return (incrementable(p.sB) && incrementable(p.sR) && 
            decrementable(p.sT) && decrementable(p.sL))
end

"""Returns true if the plaquette operator Udag can be applied to the plaquette, i.e. if the
plaquette spins are incrementable / decrementable in a loop."""
function Udag_applicable(p::Plaquette)::Bool
    return (decrementable(p.sB) && decrementable(p.sR) && 
            incrementable(p.sT) && incrementable(p.sL))
end

"""Returns the prefactor of the plaquette operator U."""
function U_prefactor(p::Plaquette)::Float64
    if U_applicable(p)
        return (increment_prefactor(p.sB) * increment_prefactor(p.sR) *
                decrement_prefactor(p.sT) * decrement_prefactor(p.sL))
    else 
        return 0.0
    end
end

"""Returns the prefactor of the plaquette operator Udag."""
function Udag_prefactor(p::Plaquette)::Float64
    if Udag_applicable(p)
        return (decrement_prefactor(p.sB) * decrement_prefactor(p.sR) *
                increment_prefactor(p.sT) * increment_prefactor(p.sL))
    else 
        return 0.0
    end
end

"""Applies the operator U to the plaquette. It returns an error if the plaquette is not
incrementable / decrementable in a loop."""
function U(p::Plaquette)
    if U_applicable(p)
        increment(p.sB)
        increment(p.sR)
        decrement(p.sT)
        decrement(p.sL)
    else
        error("U cannot be applied to the plaquette.")
    end
end

"""Applies the operator Udag to the plaquette. It returns an error if the plaquette is not
incrementable / decrementable in a loop."""
function Udag(p::Plaquette)
    if Udag_applicable(p)
        decrement(p.sB)
        decrement(p.sR)
        increment(p.sT)
        increment(p.sL)
    else
        error("Udag cannot be applied to the plaquette.")
    end
end

"""Given a plaquette, it flips it with respect the x-axis, so it exchanges the top spinZ
with the bottom SpinZ, and inverts the rungs SpinZs."""
function x_reflect(p::Plaquette)
    t = p.sT
    b = p.sB
    p.sT = b
    p.sB = t
    invert(p.sR)
    invert(p.sL)
end

"""Given a plaquette, it returns the reflected copy of the plaquette. The inversion is
with respect the x-axis, so it exchanges the top spinZ with the bottom SpinZ, and inverts
the rungs SpinZs."""
function x_reflected(p::Plaquette)
    r = copy(p)
    x_reflect(r)
    return r
end

"""Given a plaquette, it flips it with respect the y-axis, so it exchanges the right spinZ
with the left SpinZ, and inverts the top and bottom SpinZs."""
function y_reflect(p::Plaquette)
    r = p.sR
    l = p.sL
    p.sR = l
    p.sL = r
    invert(p.sT)
    invert(p.sB)
end

"""Given a plaquette, it returns the reflected copy of the plaquette. The inversion is
with respect the y-axis, so it exchanges the right spinZ with the left SpinZ, and inverts
the top and bottom SpinZs."""
function y_reflected(p::Plaquette)
    r = copy(p)
    y_reflect(r)
    return r
end
