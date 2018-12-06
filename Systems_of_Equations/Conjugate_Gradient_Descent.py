"""
This script uses conjugate gradient descent to solve a system of linear equations
conjugate gradient descent is an iterative approach which aims to minimise the error
between an initial guess and the solution in a single direction

On each iteration a new direction is chosen to minimise the error, this direction is orthogonal to all
directions used before. In this regard a solution is reached after n iterations where n is the number of equations
(Since n orthogonal vectors span R^n)
"""

import numpy as np
from matplotlib import pyplot as plt




equation_1 = lambda x: 3 - x                        # 2x + 2y = 6 -> y = 3 - x
equation_2 = lambda x: (3 - (2 * x)) / 5            # 2x + 5y = 3 -> y = (3 - 2x) / 5


def plot_gradient_descent(funcs, estimations, x_range=(-6, 6)):
    """ Uses matplotlib to plot the system of linear equations and gradient descent estimates

    func: (list[function]), List of functions to be plot over x_range
    estimations: (list[tuple]), Estimates of solution based on gradient descent
    range: (tuple) (a, b) where a is minimum x value and b in maximum x value
    """

    x_values = np.linspace(x_range[0], x_range[1], 100)        # create 100 evenly spaced points on the range x_range
    y_values_1 = [funcs[0](x) for x in x_values]               # compute equation 1 for each value of x
    y_values_2 = [funcs[1](x) for x in x_values]               # compute equation 1 for each value of x

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)

    # center the axis at (0, 0)
    ax.spines['left'].set_position('center')
    ax.spines['bottom'].set_position('center')

    estimates_x = [x[0] for x in estimations]
    estimates_y = [y[1] for y in estimations]

    plt.plot(x_values, y_values_1, color='blue', linewidth=2)         # plot equation 1
    plt.plot(x_values, y_values_2, color='red', linewidth=2)          # plot equation 1
    plt.plot(estimates_x, estimates_y, 'bo')                          # plot estimates

    # add lables, legend and title
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend(['equation 1', 'equation 2', 'estimates'], loc=2)
    plt.title("Conjugate Gradient Descent")

    # show the plot
    plt.show()


def conjugate_gradient_descent(left_hand_side, right_hand_side, initial_guess):
    """
    Uses the method of conjugate gradient descent to solve a system of linear equations
    given an initial guess for the solution.

    Solution will converge within n iterations
    """

    if len(left_hand_side) != len(right_hand_side):
        raise ValueError("dimensionality of left hand side and right hand side are not eqivilent")

    if len(left_hand_side) != len(initial_guess):
        raise ValueError("dimensionality of left hand side and the initial guess are not eqivilent")

    estimates = [initial_guess]
    current_estimate = np.copy(initial_guess)
    residual_vector = right_hand_side - left_hand_side.dot(initial_guess)
    direction_vector = np.copy(residual_vector)

    for _ in range(len(left_hand_side)):
        alpha = float(np.transpose(residual_vector).dot(residual_vector)) / float((np.transpose(direction_vector).dot(left_hand_side)).dot(residual_vector))
        current_estimate = current_estimate + alpha * direction_vector

        old_residual_vector = np.copy(residual_vector)
        residual_vector = residual_vector - alpha * (left_hand_side.dot(direction_vector))

        beta = np.transpose(residual_vector).dot(residual_vector) / np.transpose(old_residual_vector).dot(old_residual_vector)
        direction_vector = residual_vector + beta * direction_vector

        estimates.append(current_estimate)

    return estimates

def main():

    left_hand_side = np.array([[2, 2], [2, 5]])
    right_hand_side = np.array([6, 3])
    initial_guess = np.array([0, 0])

    estimates = conjugate_gradient_descent(left_hand_side, right_hand_side, initial_guess)
    plot_gradient_descent([equation_1, equation_2], estimates)

    print estimates[-1]

if __name__ == "__main__":
    main()
