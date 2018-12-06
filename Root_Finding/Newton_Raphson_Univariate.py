
"""

This script implements te Newton Raphson method for finding the roots of an equation
This is a refined version of the fixed point iteration method where the slope of g(x) is designed to be zero
This ensures high speed convergence

The algorithm works for any function f(x) = 0
an initial guess is made for the root and the tangent of the function is computed at this guess
the tangent line is followed until it bisects the x-axis
Where the tangent line bi-sects the x-axis is the next approximation of the root
This process is repeated iterativly

mathmatically this can be written as:
x_(n+1) = x_n - (f(x_n)/f'(x_n))
where f'(x_n) is the derivative

Note that the only time the Newton mehtod won't converge is when the slope is zero (i.e. at an inflection point)
In this case tangent never bisects the x-axis (we will be dividing by zero)
we can see this with the function 5 - x^2 with an initial guess of x = 0

"""

import unittest

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation


class UnitTestCases(unittest.TestCase):
    """ Class for running unit tests on functions in the Newton Raphson script

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

    def compute_function_derivative_01(self):
        """ unit tests for compute_function_02 """

        error_message = "Unit test Compute function derivative 01 Failed for input %s with solution %s"

        # check that these inputs return the correct output
        for input_val, output_val in zip([0, 0.43, -0.43, 0.1, 20], [-5.0, 0.0, 0.0, -4.73, 10795]):
            value = compute_function_derivative_01(input_val)
            self.assertTrue(np.isclose(value, output_val, rtol=0, atol=0.02, equal_nan=True), error_message % (input_val, value))

        # check that we handle inputs that are not a float or an int
        self.assertRaises(ValueError, compute_function_derivative_01, "input is a string")
        self.assertRaises(ValueError, compute_function_derivative_01, '5')

    # def test_compute_functon_02(self):


def compute_function_01(x):
    """ Computes the function 9x^3 - 5x + 20 """

    if not isinstance(x, float) and not isinstance(x, int):
        raise ValueError("Input %s does not have type float or int" % x)

    return 9*np.power(float(x), 3) - 5*float(x) + 20


def compute_function_derivative_01(x):
    """ Computes the function 27x^2 - 5 which is the derivative of 9x^3 - 5x + 20 """

    if not isinstance(x, float) and not isinstance(x, int):
        raise ValueError("Input %s does not have type float or int" % x)

    return 27*np.power(float(x), 2) - 5


def plot_function(func, x_iterations, x_range=(-6, 6)):
    """ Uses matplotlib to plot a function on a given range

    func: (function) functionn g(x) to be ploted
    range: (tuple) (a, b) where a is minimum x value and b in maximum x value
    """

    x_values = np.linspace(x_range[0], x_range[1], 100)        # create 100 evenly spaced points on the range x_range
    y_values = [func(x) for x in x_values]                     # compute g(x) for each value of x

    y_iterations = [func(x) for x in x_iterations]

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)

    # center the axis at (0, 0)
    ax.spines['left'].set_position('center')
    ax.spines['bottom'].set_position('center')

    plt.plot(x_values, y_values, color='red', linewidth=2)          # plot g(x)
    plt.plot(x_iterations, y_iterations, 'bo')

    # add lables, legend and title
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend(['g(x)', "root estimates"], loc=2)
    plt.title("Newton Raphson Method")

    # show the plot
    plt.show()


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
    plt.title("Newton Raphson Method")

    # show the plot
    plt.show()


def run_newton_raphson(initial_guess, func, derivative_func):
    """ Uses the newton raphson method to find the roots of a function f(x) given an initial guess

    initial_guess: (int/float) initial estimate of the root
    func: (function) function which computes f(x)
    derivative_func: (function) function which computes f'(x)
    """

    current_estimate = initial_guess    # initial guess for root of equation
    max_iteration = 1000                # maximum number of iterations before we give up on convergence
    num_iterations = 0                  # count the number of completed iterations
    epsilon = 0.000005                  # condition for convergence i.e. when |x_(n+1) - x_n| < epsilon

    previous_estimate = 0.0             # stores previous estimate for computing while loop (DONT INITIALISE TO SAME AS CURRENT ESTIMATE)
    iterations = [current_estimate]     # stores the root estimate for each iteration

    while np.abs(current_estimate - previous_estimate) >= epsilon:
        if compute_function_derivative_01(current_estimate) == 0.0:
            raise ValueError("Current estimate %s is an inflection point, slope is zero" % current_estimate)

        previous_estimate = current_estimate
        current_estimate = current_estimate - (compute_function_01(current_estimate) / compute_function_derivative_01(current_estimate))
        num_iterations += 1
        iterations.append(current_estimate)

        if num_iterations >= max_iteration:
            raise ValueError("Function did not converge in time")

    print "Estimate of root is: %s" % current_estimate
    print "Number of iterations is: %s" % num_iterations

    x_range = (min(min(iterations), -6), max(max(iterations), 6))
    x_range = (-20, 20)
    plot_animate_function(func, iterations, x_range)
    return current_estimate


def main():
    # unittest.main()
    run_newton_raphson(2.0, compute_function_01, compute_function_derivative_01)


if __name__ == "__main__":
    main()
