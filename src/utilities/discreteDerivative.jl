"""The following function takes an array and a separation, and computes the centered
discrete first derivative from the array"""
function discreteCenteredFirstDerivative(array, separation)
    derivative = zeros(length(array))
    len = length(array)
    for i = 1:len
        if i == 1
            derivative[i] = (array[i + 1] - array[i]) / separation
        elseif i == length(array)
            derivative[i] = (array[i] - array[i - 1]) / separation
        else
            derivative[i] = (array[i + 1] - array[i - 1]) / (2 * separation)
        end
    end
    return derivative
end

"""The following function takes an array and a separation, and computes the discrete
second derivative applying two times the discrete first derivative function."""
function discreteCenteredSecondDerivative(array, separation)
    derivative = discreteFirstDerivative(array, separation)
    derivative = discreteFirstDerivative(derivative, separation)
    return derivative
end

"""The following function takes an array and a separation, and computes the forward
discrete first derivative from the array. The resulting array will have one less element
than the original array."""
function discreteForwardFirstDerivative(array, separation)
    derivative = zeros(length(array) - 1)
    for i = 1:length(array) - 1
        derivative[i] = (array[i + 1] - array[i]) / separation
    end
    return derivative
end

"""The following function takes an array and a separation, and computes the discrete
second derivative applying two times the discrete first derivative function. The resulting
array will have two less elements than the original array."""
function discreteForwardSecondDerivative(array, separation)
    derivative = discreteForwardFirstDerivative(array, separation)
    derivative = discreteForwardFirstDerivative(derivative, separation)
    return derivative
end

"""The following function takes an array and eliminates the last element of the array."""
function eliminateLastElement(array)
    return array[1:length(array) - 1]
end

"""The following function takes an array and eliminates the first and the last element of
the array."""
function eliminateFirstAndLastElement(array)
    return array[2:length(array) - 1]
end