"""
This script uses the newton raphson method to solve a system of non-linear equations
Like with the univariate approach this is an itterative method which relies on the tangent of an n-dimensional surface
This tangent is described by the jacobian of the system of equations

The equations for Newtons Multivariate method are:
    J(x_k) * s = - F(x_k)
    x_(k+1) = x_k + s

x is the current guess for the solution
J is the jacobian of the system of linear equations
F(x) is the solution to the system of equations using the current guess

We solve for s using Guassian elimiation and then use this to update the guess for x
"""

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation
from LU_Factorisation import lu_factorisation, solve_system_of_equations


def plot_functions(estimates):
    """ uses matplotlib to plot the two non-linear functions
    also shows the locations of each estimate for the newton raphson approach

    estimates: (list of 1D vectors) List of all the estimates according to newton raphson
    """

    limit = 10

    def animate(i):
        graph.set_data(x_iterations[:i+1], y_iterations[:i+1])
        return graph

    x_range_1 = (-1.0, 1.0)                                # range for the function f1
    x_range_2 = (-limit, limit)                            # range for the function f2

    x_values_1 = np.linspace(x_range_1[0], x_range_1[1], 500)        # create 100 evenly spaced points on the range x_range
    x_values_2 = np.linspace(x_range_2[0], x_range_2[1], 500)        # create 100 evenly spaced points on the range x_range

    f1_values_a = np.sqrt(1 - np.power(x_values_1, 2))           # compute f1(x) for each value of x (top half due to square root)
    f1_values_b = -np.sqrt(1 - np.power(x_values_1, 2))          # compute f1(x) for each value of x (bottom half due to square root)
    f2_values = np.power(x_values_2, 3)                          # compute f2(x) for each value of x

    x_iterations, y_iterations = zip(*estimates)

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)

    # center the axis at (0, 0)
    ax.spines['left'].set_position('center')
    ax.spines['bottom'].set_position('center')

    # set limits on the axis
    plt.ylim((-limit, limit))   # set the ylim to bottom, top

    plt.plot(np.concatenate((x_values_1, x_values_1), axis=0), np.concatenate((f1_values_a, f1_values_b), axis=0), color='blue', linewidth=2)          # plot f1(x)
    plt.plot(x_values_2, f2_values, color='blue', linewidth=2)                                                                                         # plot f2(x)
    graph, = plt.plot([], [], 'bo', color='red')
    ani = FuncAnimation(fig, animate, frames=len(x_iterations), interval=200)

    # add lables, legend and title
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend(['f1(x)', "f2(x)"], loc=2)
    plt.title("Multivariate Newton Raphson Method")

    # show the plot
    plt.show()


def compute_jacobian(current_estimate):
    """ Computes and returns the jacobian matrix for the input values u and v
    current_estimate: (1D numpy array) the current estimate for the values of u and v as a vector

    Jacobian is computed for the equations:
        v - u^2 = 0
        u^2 + v^2 - 1 = 0
    """
    u, v = current_estimate[0], current_estimate[1]

    du_f1 = -3 * (u ** 2)               # partial derivative of first function w.r.t u
    dv_f1 = 1                           # partial derivative of first function w.r.t v

    du_f2 = 2 * u                       # partial derivative of second function w.r.t u
    dv_f2 = 2 * v                       # partial derivative of second function w.r.t v

    return np.array([[du_f1, dv_f1], [du_f2, dv_f2]], dtype=np.float)


def compute_system(current_estimate):
    """ Computes the result of the systems of equations given the guess for u and v
    current_estimate: (1D numpy array) the current estimate for the values of u and v as a vector

    Computes the equations:
        v - u^2 = 0
        u^2 + v^2 - 1 = 0
    """
    u, v = current_estimate[0], current_estimate[1]

    f1 = v - (u ** 3)
    f2 = (u ** 2) + (v ** 2) -1

    return np.array([f1, f2], dtype=np.float)


def solve_guassian_elimination(left_hand_side, right_hand_side):
    """ Use LU factorisation to solve a system of linear equations
    This is required to compute the s-vector needed for the multivariate newton raphson method
    This acts as an alternative to invering the Jacobian matrix which is computationally exspensive

    left_hand_side: (2D numpy array) The Jacobian matrix solve with the current guess for u and v
    right_hand_side: (1D numpy array) The result of the system of equations solved for the current guess of u and v
    """

    p_matrix, l_matrix, u_matrix = lu_factorisation(left_hand_side)
    return solve_system_of_equations(u_matrix, l_matrix, p_matrix.dot(right_hand_side))


def run_multivariate_newton_raphson(starting_guess):
    """ Uses the multivariate newton raphson method to solve a system of non linear equations
    starting_guess: (1D numpy array) The initial guess for the newton raphson iterative approach
    """
    current_guess = starting_guess
    estimates = [current_guess]
    for _ in range(10):
        left_hand_side = compute_jacobian(current_guess)
        right_hand_side = compute_system(current_guess)

        s_vector = solve_guassian_elimination(left_hand_side, -right_hand_side)
        current_guess = current_guess + s_vector

        estimates.append(current_guess)

    return estimates

def main():
    """

    """

    starting_guess = np.array([4, 4])

    estimates = run_multivariate_newton_raphson(starting_guess)
    plot_functions(estimates)

if __name__ == "__main__":
    main()
