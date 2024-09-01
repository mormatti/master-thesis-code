"""This module defines the SpinZ type, which represents a spin in the z axis. It is defined
by the spin value and the spin projection value in the z axis. The spin value is a positive
integer, and the spin projection value is an integer between -spin and spin. The spin value
is stored doubled, so that the spin value is 2s, where s is the actual spin value."""

mutable struct SpinZ
    ds::Int # The doubled spin value
    dsz::Int # The doubled spin projection value in the z axis

    # The constructor for SpinZ via integer (doubled spin) notation
    function SpinZ(ds::Int, dsz::Int)
        if ds < 0
            error("The spin value cannot be negative")
        end
        if ds < dsz
            error("The spin projection value cannot be larger than the spin value")
        end
        if dsz < -ds
            error("The spin projection value cannot be smaller than the negative of the spin value")
        end
        if (ds - dsz) % 2 != 0
            error("The difference between the spin value and the spin projection value must be an integer")
        end
        new(ds, dsz)
    end
end

"""The constructor for SpinZ without arguments. It returns the spin 1/2 with z-projection 1/2"""
SpinZ() = SpinZ(1, 1)

"""The constructor for SpinZ via string notation"""
SpinZ(string_ds::String, string_dsz::String) = SpinZ(ds(string_ds), ds(string_dsz))

"""Given a SpinZ object, return a deep-copy of it."""
function Base.deepcopy(spin::SpinZ)::SpinZ
    return SpinZ(spin.ds, spin.dsz)
end

"""Given an integer, it checks if it is a good double spin value"""
good_ds(ds::Int)::Bool = (ds > 0 ? true : false)

"""Given two integers, namely the doubled spin value and the doubled spin projection value,
return true if they are consistent as doubled spins. Otherwise, return false."""
good_doubled(ds::Int, dsz::Int)::Bool = (good_ds(ds) && ds >= dsz && dsz >= -ds && (ds - dsz) % 2 == 0 ? true : false)

"""Given a new value of the doubled spin projection value, it checks if it is consistent
with the doubled spin value. If it is, it assign the new value to the doubled spin projection
value. Otherwise, it throws an error."""
function set_dsz(spin::SpinZ, dsz::Int)
    if good_doubled(spin.ds, dsz)
        spin.dsz = dsz
    else
        error("The spin projection value is not consistent with the spin value")
    end
end

"""Given an integer representing the doubled spin value, return the corresponding
string notation of the spin"""
function string_spin_from_ds(ds::Int)::String
    if ds % 2 == 0
        return string(ds ÷ 2)
    else
        return string(ds, "/2")
    end
end

"""The equality operator for SpinZ. It checks if the spin and spin projection values are equal
for both SpinZ objects."""
function Base.:(==)(s1::SpinZ, s2::SpinZ)::Bool
    s1.ds == s2.ds && s1.dsz == s2.dsz
end

"""Return true if the projection value can be incremented by 1, false otherwise"""
incrementable(spin::SpinZ) = (spin.dsz < spin.ds ? true : false)

"""Return true if the projection value can be decremented by 1, false otherwise"""
decrementable(spin::SpinZ) = (spin.dsz > -spin.ds ? true : false)

"""Increment the projection value by 1. If the projection value is already the maximum value,
an error is thrown."""
function increment(spin::SpinZ)
    if incrementable(spin)
        spin.dsz += 2
    else
        error("The spin projection value cannot be incremented")
    end
end

"""Return the matrix element of the spin ladder operator S+ between the spin and its
incremented value"""
function increment_prefactor(spin::SpinZ)
    if incrementable(spin)
        s = spin.ds / 2
        m = spin.dsz / 2
        return sqrt(s * (s + 1) - m * (m + 1))
    else
        return 0
    end
end

"""Return the matrix element of the spin ladder operator S- between the spin and its
decremented value"""
function decrement_prefactor(spin::SpinZ)
    if decrementable(spin)
        s = spin.ds / 2
        m = spin.dsz / 2
        return sqrt(s * (s + 1) - m * (m - 1))
    else
        return 0
    end
end

"""Decrement the projection value by 1. If the projection value is already the minimum value,
an error is thrown."""
function decrement(spin::SpinZ)
    if decrementable(spin)
        spin.dsz -= 2
    else
        error("The spin projection value cannot be decremented")
    end
end

"""This function converts a spin-s string notation of spin to its doubled value ds"""
function ds(string::String)::Int
    if length(string) == 1
        return 2 * parse(Int, string)
    elseif length(string) == 3
        return parse(Int, string[1])
    elseif length(string) == 4
        return -parse(Int, string[2])
    else
        error("Invalid spin notation")
    end
end

"""This function converts a sz spin string notation of spin to its doubled value dsz"""
dsz(string::String)::Int = ds(string)

"""Given a doubled spin value, returns a SpinZ with a random spin projection value. The
spin projection value is chosen uniformly at random from the set of possible values."""
function SpinZ_random(ds::Int)::SpinZ
    if good_ds(ds)
        return SpinZ(ds, rand(-ds:2:ds))
    else
        error("The spin value is not valid")
    end
end

"""Given a doubled spin value, return the list of all possible SpinZ objects with that doubled
spin value."""
function SpinZ_list(ds::Int)::Vector{SpinZ}
    if good_ds(ds)
        return [SpinZ(ds, dsz) for dsz in -ds:2:ds]
    else
        error("The spin value is not valid")
    end
end

"""It flips the value of the projection of the spin."""
function invert(spin::SpinZ)
    spin.dsz = -spin.dsz
end