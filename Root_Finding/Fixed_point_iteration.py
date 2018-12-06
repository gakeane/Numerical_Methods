
"""
a fixed point of a function is an element in the functions domain that is mapped to itself
i.e. c is a fixed point of a function f(x) if f(c) = c
By converting a function f(x) = 0 to the form g(x) = x we can use the fixed point iteration algorithm
to find the roots of g(x) = x which in turn will be the roots of f(x) = 0


The fixed point iteration algorithm is as follows

Given a function f(x) = 0
Convert it to the form g(x) = x
Choose an initial guess for the fixed point x

Note the fixed point iteration algorithm will only converge if |g'(r)| < 1
Where g' is the derivate of g and r is a root of g
i.e. we will only see convergence for roots where the slope of the function is less than 1
The closer the slope is to zero at the root the faster the algorithm will converge

In this script we give an example of using the fixed point algorithm to find the roots of the function x^4 - x - 10 = 0
"""

import unittest

import numpy as np
from matplotlib import pyplot as plt


class UnitTestCases(unittest.TestCase):
    """ Class for running unit tests on functions in the fixed point itteration script

    The names of any unit tests should be prefixed with "test"
    The unit test frame work consist of 4 checks assertEqual, asertTrue, assertFalse, assertRaises
    we run all the unit tests by calling unittest.main()
    """

    def test_run_fixed_point_iteration(self):
        """ unit test for run_fixed_point_iteration """

        error_message = "Unit test fixed point iteration function Failed for input %s with solution %s"

        root = run_fixed_point_iteration(2.0, compute_function_01)
        self.assertTrue(np.isclose(root, 1.85558, rtol=0, atol=0.002, equal_nan=True), error_message % (2.0, root))

        root = run_fixed_point_iteration(2.0, compute_function_02)
        self.assertTrue(np.isclose(root, 0.6437, rtol=0, atol=0.002, equal_nan=True), error_message % (2.0, root))

    def test_compute_functon_01(self):
        """ unit tests for compute_function_01 """

        error_message = "Unit test Compute function 01 Failed for input %s with solution %s"

        # check that these inputs return the correct output
        for input_val, output_val in zip([0, 2.3, 1.5, 300, -5, -10], [1.77827941, 1.872734787, 1.8415116, 4.1960477, 1.4953487, 0.0]):
            self.assertTrue(np.isclose(compute_function_01(input_val), output_val, rtol=0, atol=0.00002, equal_nan=True), error_message % (input_val, compute_function_01(input_val)))

        # check that values less than -10 throw an exception since they will result in complex numbers
        self.assertRaises(ValueError, compute_function_01, -13)
        self.assertRaises(ValueError, compute_function_01, -12.0)

        self.assertRaises(ValueError, compute_function_01, "input is a string")
        self.assertRaises(ValueError, compute_function_01, '5')


def compute_function_01(x):
    """ Computes the function (x + 10)^(1/4) converted from f(x) = x^4 - x - 10"""

    if not isinstance(x, float) and not isinstance(x, int):
        raise ValueError("Input %s does not have type float or int" % x)

    if x < -10:
        raise ValueError("Result is a complex number for input %s" % x)

    return np.power((float(x) + 10.0), (0.25))


def compute_function_02(x):
    """ Computes the function cos(x)^2 , converted from f(x) = ((cos(x)^2)/x) - 1"""

    if not isinstance(x, float) and not isinstance(x, int):
        raise ValueError("Input %s does not have type float or int" % x)

    return np.power(np.cos(x), 2)


def plot_function(func, iterations_x, x_range=(-6, 6)):
    """ Uses matplotlib to plot a function on a given range

    func: (function) functionn g(x) to be ploted
    range: (tuple) (a, b) where a is minimum x value and b in maximum x value
    """

    x_values = np.linspace(x_range[0], x_range[1], 100)        # create 100 evenly spaced points on the range x_range
    y_values = [func(x) for x in x_values]                     # compute g(x) for each value of x

    y_iterations = [func(x) for x in iterations_x]

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)

    # center the axis at (0, 0)
    ax.spines['left'].set_position('center')
    ax.spines['bottom'].set_position('center')

    plt.plot(x_values, y_values, color='red', linewidth=2)          # plot g(x)
    plt.plot(x_values, x_values, color='blue', linewidth=1)         # plot x = y
    plt.plot(iterations_x, y_iterations, 'bo')

    # add lables, legend and title
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend(['g(x)', 'x = y'], loc=2)
    plt.title("Fixed point iteration")

    # show the plot
    plt.show()


def run_fixed_point_iteration(initial_guess, func):
    """
    Runs the fixed point algorithm on a function g(x) with an initial guess
    convergence condition is set at 0.000005
    will allow up to 200 iterations

    initial guess: (int/float) the initial guess for the root
    func: (function) function which calculates g(x)
    """

    current_estimate = initial_guess    # initial guess for root of equation
    max_iteration = 10000               # maximum number of iterations before we give up on convergence
    num_iterations = 0                  # count the number of completed iterations
    epsilon = 0.000005                  # condition for convergence i.e. when |x_(n+1) - x_n| < epsilon

    previous_estimate = 0.0             # stores previous estimate for computing while loop (DONT INITIALISE TO SAME AS CURRENT ESTIMATE)
    iterations = [current_estimate]     # stores the root estimate for each iteration

    while np.abs(current_estimate - previous_estimate) >= epsilon:
        previous_estimate = current_estimate
        current_estimate = func(current_estimate)
        num_iterations += 1
        iterations.append(current_estimate)

        if num_iterations >= max_iteration:
            raise ValueError("Function did not converge in time")

    print "Estimate of root is: %s" % current_estimate
    print "Number of iterations is: %s" % num_iterations

    plot_function(func, iterations)
    return current_estimate


def main():
    unittest.main()           # uncomment to run unit tests
    # run_fixed_point_iteration(2.0, compute_function_02)


if __name__ == "__main__":
    main()
