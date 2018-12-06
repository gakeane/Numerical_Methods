"""
This script implements the secant method for approximating the roots of a function
The secant method is similar to the Newton Raphson method except that instead of the tangent line
We use the line that intercepts the previous two guesses
Where this line intercepts the x axis is the new guess

the secant line is an approximation of the tangent line
Mathmatically this look like:

f'(x_n) ~= (f(x_n) - f(x_(n-1))) / (x_n - x_(n-1))

This is the equation for calculating the slope of a line
Subbing this into the Newton raphson equation gives

x_n = x_(n-1) - f(x_(n-1)) * (x_n - x_(n-1)) / (f(x_n) - f(x_(n-1)))

This approach removes the issue of the slope being zero
"""

import unittest

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation


class UnitTestCases(unittest.TestCase):
    """ Class for running unit tests on functions in the Secanrscript

    The names of any unit tests should be prefixed with "test"
    The unit test frame work consist of 4 checks assertEqual, asertTrue, assertFalse, assertRaises
    we run all the unit tests by calling unittest.main()
    """

    def test_compute_functon_01(self):
        """ unit tests for compute_function_01 """

        error_message = "Unit test Compute function 01 Failed for input %s with solution %s"

        # check that these inputs return the correct output
        for input_val, output_val in zip([0, 1.0, -1.5, -1.446366], [20.0, 24.0, -2.875, 0.0]):
            value = compute_function_01(input_val)
            self.assertTrue(np.isclose(value, output_val, rtol=0, atol=0.002, equal_nan=True), error_message % (input_val, value))

        # check that we handle inputs that are not a float or an int
        self.assertRaises(ValueError, compute_function_01, "input is a string")
        self.assertRaises(ValueError, compute_function_01, '5')


def compute_function_01(x):
    """ Computes the function 9x^3 - 5x + 20 """

    if not isinstance(x, float) and not isinstance(x, int):
        raise ValueError("Input %s does not have type float or int" % x)

    return 9*np.power(float(x), 3) - 5*float(x) + 20


def plot_animate_function(func, x_iterations, x_range=(-6, 6)):
    """ Uses matplotlib to plot a function on a given range with animations

    func: (function) functionn g(x) to be ploted
    range: (tuple) (a, b) where a is minimum x value and b in maximum x value
    """

    def animate(i):
        graph.set_data(x_iterations[:i+1], y_iterations[:i+1])
        return graph

    x_values = np.linspace(x_range[0], x_range[1], 100)        # create 100 evenly spaced points on the range x_range
    y_values = [func(x) for x in x_values]                     # compute g(x) for each value of x

    y_iterations = [func(x) for x in x_iterations]

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)

    # center the axis at (0, 0)
    ax.spines['left'].set_position('center')
    ax.spines['bottom'].set_position('center')

    plt.plot(x_values, y_values, color='red', linewidth=2)          # plot g(x)
    graph, = plt.plot([], [], 'bo')
    ani = FuncAnimation(fig, animate, frames=len(x_iterations), interval=200)

    # add lables, legend and title
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend(['g(x)', "root estimates"], loc=2)
    plt.title("Secant Method")

    # show the plot
    plt.show()


def run_secant_iteration(initial_guess_1, initial_guess_2, func):
    """ Uses the secant method to find the roots of a function f(x) given two initial guesses
    The first secant line will be drawn between the two initial guesses

    initial_guess_1: (int/float) first initial estimate of the root
    initial_guess_2: (int/float) second initial estimate of the root
    func: (function) function which computes f(x)
    """

    previous_estimate = initial_guess_1   # first initial guess for root of equation
    current_estimate = initial_guess_2    # second initial guess for root of equation

    max_iteration = 5000                  # maximum number of iterations before we give up on convergence
    num_iterations = 0                    # count the number of completed iterations
    epsilon = 0.000005                    # condition for convergence i.e. when |x_(n+1) - x_n| < epsilon

    iterations = [previous_estimate, current_estimate]     # stores the root estimate for each iteration#

    if np.abs(current_estimate - previous_estimate) < epsilon:
        raise ValueError("Initial Guess 1: %s and Initial Guess 2: %s are too close" % (initial_guess_1, initial_guess_2))

    while np.abs(current_estimate - previous_estimate) >= epsilon:

        secant_slope = (compute_function_01(current_estimate) - compute_function_01(previous_estimate)) / (current_estimate - previous_estimate)

        previous_estimate = current_estimate
        current_estimate = current_estimate - (compute_function_01(current_estimate) / secant_slope)
        num_iterations += 1
        iterations.append(current_estimate)

        if num_iterations >= max_iteration:
            raise ValueError("Function did not converge in time")

    print "Estimate of root is: %s" % current_estimate
    print "Number of iterations is: %s" % num_iterations

    x_range = (-20, 20)
    plot_animate_function(func, iterations, x_range)
    return current_estimate


def main():
    # unittest.main()
    run_secant_iteration(2.0, 2.2, compute_function_01)


if __name__ == "__main__":
    main()
