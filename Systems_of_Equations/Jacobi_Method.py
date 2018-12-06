"""
This script implements the Jacobi method for solving systems of linear equations
The Jacobi method is an extension of the fixed point method
Given a system of N equations with N unknowns we solve the ith equation with the ith unknown
I.e the equations 3u + v = 5 and u + 2v = 5 got to u = (5 - v) / 3 and v = (5 - u) / 2

We take an initial gues for u and v (i.e. x = [u_0, v_0])
Ax = b can be rewritten as (D + L + U)x = b
where D is the diagonal of A
U is the upper triangular matrix of A where the diagonal is zero
L is the lower traingular matrix of A where the diagonal is zero

We can then rewrite the above equations as
x_(n+1) = inv(D) * (b - (L + U) * x_n)
note that since D is a diagonal matrix the inverse is just a diagonal matrix of it's recipricals

The jacobi method will not always converge to a solution, a sufficinet condition for convergence is if the matrix of equations is
diagonal dominant (i.e. the value of the diagonal is greater than the sum of all the other values in the row)
"""

import numpy as np


def _check_matrix_shape(matrix):
    """ Confrims that a matrix is an NxN matrix """

    matrix_shape = matrix.shape

    if len(matrix_shape) != 2:
        raise ValueError("Matrix must be two dimensionoal")

    if matrix_shape[0] != matrix_shape[1]:
        raise ValueError("Matrix is not square")


def is_diagonal_dominant(matrix):
    """ returns True if an nxn matrix is diagonaly dominant, False otherwise
    A diagonally dominant matrix is one where the diagonal element of each row is greater than the sum of all other elements in the row

    matrix: (2D nummpy array) An nxn matrix
    """

    _check_matrix_shape(matrix)

    for row in range(len(matrix)):

        # check if any of the rows are not diagonal dominant
        if np.abs(matrix[row][row]) <= np.sum(np.abs(matrix[row])) - np.abs(matrix[row][row]):
            return False

    return True


def DLU_reduction(matrix):
    """ Given an nxn matrix A returns the matricies D L and U
    D is the diagonal matrix of A
    U is the upper triangular matrix of A with the diagonal set to zero
    L is the lower triangular matrix of A with the diagonal set to zero

    matrix: (2D nummpy array) An nxn matrix
    """

    _check_matrix_shape(matrix)

    N = len(matrix)
    D = np.zeros((N, N))
    U = np.zeros((N, N))
    L = np.zeros((N, N))

    for i in range(N):
        # create the diagonal matrix
        D[i][i] = matrix[i][i]

        # create upper triangular matrix
        for j in range(i + 1, N):
            U[i][j] = matrix[i][j]

        # create lower triangualr matrix
        for j in range(i):
            L[i][j] = matrix[i][j]

    return D, U, L


def invert_diagonal_matrix(matrix):
    """ Given a diagonal matrix returns the inverted matrix """

    _check_matrix_shape(matrix)
    matrix_copy = np.copy(matrix)           # make a copy to aviod clobbering the original matrix

    for i in range(len(matrix)):
        if matrix_copy[i][i] == 0:
            raise ValueError("One of the diagonal elements is zero, matrix is not invertable")

        matrix_copy[i][i] = 1 / matrix_copy[i][i]
    return matrix_copy


def run_jacobi_method(left_hand_side, right_hand_side, initial_guess):
    """ Uses the Jacobi method for determining the roots of a system of equations

    left_hand_side: (numpy 2D array) NxN matrix representing the left hand side of the system of equations
    right_hand_side: (numpy 1D array) Nx1 vector representing the right hand side of the system of equations
    initial guess: (numpy 1D array) Nx1 vector representing the initial guess for teh roots of the system of equations
    """

    _check_matrix_shape(left_hand_side)

    if len(right_hand_side.shape) != 1:
        raise ValueError("The right hand side of the system must be specified as a vector")

    if len(initial_guess.shape) != 1:
        raise ValueError("The initial guess must be specified as a vector")

    if len(left_hand_side) != len(right_hand_side):
        raise ValueError("left hand side and right hand side of system have incompatible sizes")

    if len(left_hand_side) != len(initial_guess):
        raise ValueError("system of equations and initial guess have incompatible sizes")

    current_guess = initial_guess                           # initalise as the initial guess
    previous_guess = np.ones(len(left_hand_side)) * 99      # we need to  avoid this being initalised to the same as the current guess
    epsilon = 0.000005                                      # stopping criteria when difference sum(abs(current_guess - abs(previos_guess))) < epsilon

    max_iterations = 20000
    num_iterations = 0
    iterations = [current_guess]

    D, U ,L = DLU_reduction(left_hand_side)
    D_inv = invert_diagonal_matrix(D)

    while np.sum(np.abs(previous_guess - current_guess)) > epsilon:

        if num_iterations >= max_iterations:
            raise ValueError("Jocabi Method did not converge after 20000 iterations")

        previous_guess = current_guess
        current_guess = np.dot(D_inv, (right_hand_side - np.dot((L + U), current_guess)))
        iterations.append(current_guess)
        num_iterations += 1

    print current_guess
    print num_iterations


def main():

    left_hand_side = np.array([[3, -1, 0, 0, 0, 0.5], [-1, 3, -1, 0, 0.5, 0], [0, -1, 3, -1, 0, 0], [0, 0, -1, 3, -1, 0], [0, 0.5, 0, -1, 3, -1], [0.5, 0, 0, 0, -1, 3]])
    right_hand_side = np.array([2.5, 1.5, 1, 1, 1.5, 2.5])
    initial_guess = np.array([0, 0, 0, 0, 0, 0])

    run_jacobi_method(left_hand_side, right_hand_side, initial_guess)


if __name__ == "__main__":
    main()
