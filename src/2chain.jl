"""The chain type. It has a single field, a vector of plaquettes. The plaquettes are ordered
from left to right."""

include("0spinZ.jl")
include("1plaquette.jl")


using LinearAlgebra

mutable struct Chain
    plaquettes::Vector{Plaquette}
    closed::Bool

    """The constructor for chain. It takes a vector of plaquettes."""
    Chain(plaquettes::Vector{Plaquette}, closed) = new(plaquettes, closed)
end

"""The constructor for the chain. It takes no arguments and returns the void chain.
The void chain is a chain with no plaquettes and it is closed."""
Chain() = Chain(Plaquette[], false)

"""It returns a deep-copy of the chain."""
function Base.deepcopy(chain::Chain)::Chain
    return Chain(deepcopy(chain.plaquettes), deepcopy(chain.closed))
end

"""Returns true if the chain is a void chain."""
function is_empty(chain::Chain)::Bool
    return length(chain.plaquettes) == 0
end

"""It returns the number of plaquettes in the chain."""
function Base.length(chain::Chain)::Int
    return length(chain.plaquettes)
end

"""The constructor for the chain. It takes a single plaquette and returns a chain with
that plaquette as its only element. The chain is by default not closed."""
Chain(plaquette::Plaquette) = Chain([plaquette], false)

"""It returns the first plaquette of the chain."""
first(chain::Chain)::Plaquette = chain.plaquettes[1]

"""It returns the last plaquette of the chain."""
last(chain::Chain)::Plaquette = chain.plaquettes[end]

"""The "precedes" binary operator, it returns true if the common rung between the two
chains has the same spin"""
function Base.:>(c1::Chain, c2::Chain)::Bool
    if is_empty(c1) || is_empty(c2)
        return true
    else 
        return last(c1) > first(c2)
    end
end

"""The concatenation operator for chain. It takes two chains and returns a new chain
(notice, that new is important, it returns a copy of the two chains concatenated) with 
the plaquettes of the two chains concatenated."""
function Base.:(*)(c1::Chain, c2::Chain)::Chain
    if c1.closed || c2.closed
        error("Cannot concatenate closed chains")
    else
        if is_empty(c1) && is_empty(c2)
            return Chain()
        elseif is_empty(c1)
            return deepcopy(c2)
        elseif is_empty(c2)
            return deepcopy(c1)
        else
            pA, pB = last(c1), first(c2)
            if pA > pB
                return Chain(vcat(deepcopy(c1.plaquettes), deepcopy(c2.plaquettes)), false)
            else
                error("Cannot concatenate chains with incompatible rungs")
            end
        end
    end
end

"""The equality operator for chain. It checks if the two chains have the same plaquettes
in the same order."""
Base.:(==)(c1::Chain, c2::Chain)::Bool = 
    (length(c1.plaquettes) == length(c2.plaquettes) 
    && (c1.closed == c2.closed)
    && all(c1.plaquettes .== c2.plaquettes))

"""The "gauge-precedes" binary operator, it returns true if the common rung between the two chains
has the same spin, and the nodes are consistent with the Gauss' law"""
▷(c1::Chain, c2::Chain)::Bool = 
    is_empty(c1) || is_empty(c2) || last(c1) ▷ first(c2)

"""The "weak precedes" binary operator. It returns true if the longitudinal polarization of 
the neearest neighbouring plaquettes of the two chains are the same."""
Base.:≥(c1::Chain, c2::Chain)::Bool = 
    is_empty(c1) || is_empty(c2) || last(c1) ≥ first(c2)

"""The "gauge-concatenation operator". If the two nearest neighbouring plaquettes has the
same longitudinal polarization, it adjusts the common rungs of the two chains in such a way
that the chain are concatenating and compatible with the Gauss' law. If the two plaquettes
have different longitudinal polarization, it throws an error."""
function ⊕(c1::Chain, c2::Chain)::Chain
    if c1.closed || c2.closed
        error("Cannot concatenate closed chains")
    else
        if is_empty(c1) && is_empty(c2)
            return Chain()
        elseif is_empty(c1)
            return deepcopy(c2)
        elseif is_empty(c2)
            return deepcopy(c1)
        else
            pA, pB = last(c1), first(c2)
            if pA ≥ pB
                c1p = deepcopy(c1)
                c2p = deepcopy(c2)
                dszCommon = pB.sT.dsz - pA.sT.dsz
                last(c1p).sR.dsz = dszCommon
                first(c2p).sL.dsz = dszCommon
                return Chain(vcat(c1p.plaquettes, c2p.plaquettes), false)
            else
                error("Cannot gauge-concatenate chains with different longitudinal polarization")
            end
        end
    end
end

"""Returns a random chain of a given length L."""
function random_Chain(L::Int, dsh::Int, dsv::Int)::Chain
    if L <= 0
        return Chain()
    elseif L == 1
        return Chain(Plaquette_random(dsh, dsv))
    else
        chain = Chain(Plaquette_random(dsh, dsv))
        for i in 2:L
            rp = Chain(Plaquette_random(dsh, dsv))
            set_dsz(rp.plaquettes[1].sL, last(chain).sR.dsz)
            chain = chain * rp
        end
        return chain
    end
end

"""Returns a random closed chain of a given length L."""
function random_closed_Chain(L::Int, dsh::Int, dsv::Int)::Chain
    chain = random_Chain(L, dsh, dsv)
    set_dsz(first(chain).sL, last(chain).sR.dsz)
    chain.closed = true
    return chain
end

"""Prints the chain in a readable way."""
function Base.show(chain::Chain; io::IO = Base.stdout)
    if is_empty(chain)
        print(io, "Void chain")
    else
        print("∘∘")
        for i in 1:length(chain.plaquettes)
            print(sign_string(chain.plaquettes[i].sT.dsz))
            print("∘∘")
        end
        println()
        for i in 1:length(chain.plaquettes)
            print(sign_string(chain.plaquettes[i].sL.dsz))
            print("  ")
        end
        println(sign_string(chain.plaquettes[end].sR.dsz))
        print("∘∘")
        for i in 1:length(chain.plaquettes)
            print(sign_string(chain.plaquettes[i].sB.dsz))
            print("∘∘")
        end
        println()
    end
end

"""Converts an integer in string, adding a + sign if positive, and - if
negative."""
function sign_string(n::Int)::String
    if n >= 0
        return "+" * string(n)
    else
        return string(n)
    end
end

"""After checking if all the nearest neighbours plaquettes have the common rung with same
spin, it checks if the chain is consistent with the Gauss' law and returns true if this is
the case."""
function is_chain_gauss(chain::Chain)::Bool
    if chain == Chain()
        return true
    else
        pl = chain.plaquettes
        for i in 1:length(pl)-1
            if pl[i] ▷ pl[i+1]
                return false
            end
        end
        return true
    end
end

"""It returns true if the chain can be closed, i.e. if the first and last plaquettes are
compatible for concatenation."""
function lockable(chain::Chain)::Bool
    if is_empty(chain)
        return true
    else
        return last(chain) > first(chain)
    end
end

"""It returns true if the chain can be gauge closed, i.e. if the first and last plaquettes
are compatible for gauge concatenation."""
function gauge_lockable(chain::Chain)::Bool
    if is_empty(chain)
        return true
    else
        return last(chain) ▷ first(chain)
    end
end

"""It closes the chain if possible, so if the initial and final rungs have the same spin.
It returns nothing if the chain cannot be closed."""
function lockChain(chain::Chain)
    if lockable(chain)
        chain.closed = true
    else
        error("The chain cannot be closed.")
    end
end

"""It returns true if the chain can be gauge-closed, i.e. if the first and last plaquettes are
compatible for gauge-concatenation."""
function gauge_lockable(chain::Chain)::Bool
    if is_empty(chain)
        return true
    else
        return last(chain) ≥ first(chain)
    end
end

"""It gauge-closes the chain if possible. It returns an error if the chain cannot be closed."""
function gauge_lock(chain::Chain)
    if gauge_lockable(chain)
        chain.closed = true
        if length(chain) ≥ 1
            dszCommon = first(chain).sT.dsz - last(chain).sT.dsz
            first(chain).sL.dsz = dszCommon
            last(chain).sR.dsz = dszCommon
        end
    else
        error("The chain cannot be gauge-closed.")
    end
end

"""It opens the chain. It does not anything if the chain is already opened."""
function openChain(chain::Chain)
    chain.closed = false
end

"""Returns true if all the plaquettes of the chain have the same longitudinal polarization,
which is a given number"""
function belongs_to_longpol_sector(chain::Chain, sector::Int)::Bool
    if is_empty(chain)
        return true
    else
        for p in chain.plaquettes
            if longitudinal_polarization(p) ≠ sector
                return false
            end
        end
        return true
    end
end

"""It adjusts the left rung of the chain, i.e. given a position i it substitutes the spin 
of the left rung of the plaquette at at position i+1 with the spin of the right rung of 
the plaquette at i. It does nothing if the chain is empty or if the index is 1 and the chain is
closed."""
function adjust_right_rung(chain::Chain, index::Int)
    if is_empty(chain)
        return
    elseif index == length(chain)
        if chain.closed
            set_dsz(first(chain).sL, last(chain).sR.dsz)
        end
    else
        set_dsz(chain.plaquettes[index+1].sL, chain.plaquettes[index].sR.dsz)
    end
end

"""It adjusts the right rung of the chain, i.e. given a position i it substitutes the spin
of the right rung of the plaquette at at position i-1 with the spin of the left rung of
the plaquette at i. It does nothing if the chain is empty or if the index is 1 and the chain is
closed."""
function adjust_left_rung(chain::Chain, index::Int)
    if is_empty(chain)
        return
    elseif index == 1
        if chain.closed
            set_dsz(last(chain).sR, first(chain).sL.dsz)
        end
    else
        set_dsz(chain.plaquettes[index-1].sR, chain.plaquettes[index].sL.dsz)
    end
end

"""It adjusts the rungs of the chain. It does nothing if the chain is empty or if the 
index is 1 and the chain is closed."""
function adjust_rungs(chain::Chain, index::Int)
    adjust_right_rung(chain, index)
    adjust_left_rung(chain, index)
end

function U_applicable(chain::Chain, index::Int)::Bool
    return U_applicable(chain.plaquettes[index])
end

function Udag_applicable(chain::Chain, index::Int)::Bool
    return Udag_applicable(chain.plaquettes[index])
end

"""It applies the U operator to the plaquette at position index of the chain."""
function U(chain::Chain, index::Int)
    if is_empty(chain)
        return
    elseif index == 1 && length(chain) == 1 && chain.closed
        error("The chain is a single plaquette closed. It cannot be transformed by U.")
    else
        U(chain.plaquettes[index])
        adjust_rungs(chain, index)
    end
end

"""It applies the Udag operator to the plaquette at position index of the chain."""
function Udag(chain::Chain, index::Int)
    if is_empty(chain)
        return
    elseif index == 1 && length(chain) == 1 && chain.closed
        error("The chain is a single plaquette closed. It cannot be transformed by Udag.")
    else
        Udag(chain.plaquettes[index])
        adjust_rungs(chain, index)
    end
end

"""The prefactor of the U operator applied to the plaquette at position index of the chain."""
U_prefactor(chain::Chain, index::Int) = U_prefactor(chain.plaquettes[index])

"""The prefactor of the Udag operator applied to the plaquette at position index of the chain."""
Udag_prefactor(chain::Chain, index::Int) = Udag_prefactor(chain.plaquettes[index])

"""It translates the chain to left by 1 unit (1 plaquette). If the chain is not closed,
it returns an error. It takes an optional index which represents the number of plaquettes
to translate. If the index is not given, it is set to 1. If the index is negative, it means
a translation to the right. If the index is greater than the length of the chain, it is computed
modulo the length of the chain. It uses the function circshift to perform the translation."""
function translate(chain::Chain, index::Int=1)
    if !chain.closed
        error("The chain is not closed. It cannot be translated.")
    else
        circshift!(chain.plaquettes, index)
    end
end

"""It returs the translated chain to left using the function translate."""
function translated(chain::Chain, index::Int=1)
    newChain = deepcopy(chain)
    translate(newChain, index)
    return newChain
end

"""It reflects the chain with respect to the x axis. It does nothing if the chain is empty."""
function x_reflect(chain::Chain)
    for p in chain.plaquettes
        x_reflect(p)
    end
end

"""It returns the reflected chain with respect to the x axis using the function x_reflect."""
function x_reflected(chain::Chain)
    newChain = deepcopy(chain)
    x_reflect(newChain)
    return newChain
end

"""It reflects the chain with respect to the y axis. It does nothing if the chain is empty."""
function y_reflect(chain::Chain)
    reverse!(chain.plaquettes)
    for p in chain.plaquettes
        y_reflect(p)
    end
end

"""It returns the reflected chain with respect to the y axis using the function y_reflect."""
function y_reflected(chain::Chain)
    newChain = deepcopy(chain)
    y_reflect(newChain)
    return newChain
end