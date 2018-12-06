"""
This script implements the bi-scetion method for finding roots
This method is slow but is certain to find a root for a continous function if one exists

given an interval [a, b] such that f(a) < 0 and f(b) > 0 there is certain to be at least one root on the interval
We divide this interval in half such that c = (a + b)/2
if f(r) = 0 we found the root
if f(c) < 0 we assign a = c and b stays the same
if f(c) > 0 we assign b = c and a stays the same
We continue dividing the interval until abs(a - b) < e
The approximation of the root is than (a + b) / 2
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


def run_bisection_iteration(a, b, func):
    """ Uses the secant method to find the roots of a function f(x) given two initial guesses
    The first secant line will be drawn between the two initial guesses

    a: (int/float) low value for the interval
    b: (int/float) high value for the interval
    func: (function) function which computes f(x)
    """

    if a >= b:
        raise ValueError("low interval range %s must be less than high interval range %s" % (a, b))

    if np.sign(compute_function_01(a)) == np.sign(compute_function_01(b)):
        if np.sign(compute_function_01(a)) == 1:
            raise ValueError("Function is positive at both ends of the interval")
        if np.sign(compute_function_01(a)) == -1:
            raise ValueError("Function is negative at both ends of the interval")


    low_interval = a            # first initial guess for root of equation
    high_interval = b           # second initial guess for root of equation

    max_iteration = 5000        # maximum number of iterations before we give up on convergence
    num_iterations = 0          # count the number of completed iterations
    epsilon = 0.000005          # condition for convergence i.e. when |x_(n+1) - x_n| < epsilon

    iterations = [(a + b) / 2]     # stores the root estimate for each iteration#

    while np.abs(b - a) >= epsilon:

        c = (a + b) / 2

        # special (unlikely) case where we find the root exactly
        if np.sign(compute_function_01(c)) == 0:
            a = c
            b = c

        # if f(c) has the same sign as f(a) then a = c
        if np.sign(compute_function_01(c)) == np.sign(compute_function_01(a)):
            a = c

        # if f(c) has the same sign as f(b) then b = c
        if np.sign(compute_function_01(c)) == np.sign(compute_function_01(b)):
            b = c

        num_iterations += 1
        iterations.append(c)

        if num_iterations >= max_iteration:
            raise ValueError("Function did not converge in time")

    print "Estimate of root is: %s" % c
    print "Number of iterations is: %s" % num_iterations

    x_range = (-20, 20)
    plot_animate_function(func, iterations, x_range)
    return c


def main():
    # unittest.main()
    run_bisection_iteration(-20.0, 20.0, compute_function_01)


if __name__ == "__main__":
    main()
