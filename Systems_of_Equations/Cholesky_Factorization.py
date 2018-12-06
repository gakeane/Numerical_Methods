"""
This script use cholesky factorisation to solve a system of linear equations where the left hand side
of the system of equations from a symmetric positive definite matrix
A matrix is symmetric if transpose(A) = A
A matrix is positive definite if transpose(x)Ax > 0 for all vecotes x != 0

Note that a positive definite matrix will always be non-singular (i.e. invertable)
A symmetric matrix will be positive definite if and only if all it's eigen values are positive
The determiniant of a matrix is the product of it's eigen values so will always be positive for a symmetric positive definite matrix
"""

import numpy as np


def cholesky_factorization(matrix):
    """ Returns the upper triangular matrix R resulting from a cholesky factorisation on the input matrix A
    The cholesky factorization is such that a = transpose(R)R
    The cholesky_factorization only works on matricies that are symmetric positive definite
    Cholesky factorization is done recursively on the (n-1) x (n-1) lower right matrix until n = 1

    matrix: (numpy 2D array) symmetric positive definite nxn matrix
    """

    matrix_copy = np.copy(matrix)               # make a copy of the matrix to prevent clobbering the existing matrix

    for i in range(len(matrix_copy)):

        # calculate the top left corner of the sub-matrix
        matrix_copy[i][i] = np.sqrt(matrix_copy[i][i])

        # calculate the top right side of the sub-matrix
        matrix_copy[i][i + 1:] = matrix_copy[i][i + 1:] / matrix_copy[i][i]

        # calculate the bottom left side of the sub-matrix
        for row in range(i + 1, len(matrix_copy)):
            matrix_copy[row][i] = 0

        # subtract outer product of transpose(u)u from sub matrix
        matrix_copy[i + 1:, i + 1:] = matrix_copy[i + 1:, i + 1:] - np.outer(matrix_copy[i][i + 1:], matrix_copy[i][i + 1:])

    return matrix_copy, np.transpose(matrix_copy)


def solve_system_of_equations(upper_triangular, lower_triangular, right_hand_side):
    """ Solves a system of linear equations using back substitution where the matrix A has been
    factorised into an upper triangular matrix and a lower triangular matrix
    We first solve Lc = b and the Ux = c where c, b and x are nx1 vectors and x = right_hand_side of system

    upper_triangular: (numpy 2D array), upper triangular nxn matrix
    lower_triangular: (numpy 2D array), lower triangular nxn matrix
    right_hand_side: (numpy 2D array), nx1 vector representing the right hand side of the system of equations
    """

    # use back substituion with the lower triangular matrix and the right hand side to solve for the c_vector
    c_vector = [1] * len(lower_triangular)
    for row in range(len(lower_triangular)):

        # substitute in the known c_vector values for all the non-zero coefficients and subtract from the right hand side for the current row
        for c_index in range(row):
            right_hand_side[row] = right_hand_side[row] - c_vector[c_index] * lower_triangular[row][c_index]

        # calculate the c_vector value for this row with the updated right hand side value
        c_vector[row] = right_hand_side[row] / lower_triangular[row][row]


    # use back substitution with the upper triangular matrix and the c_vector to solve for the roots
    roots = [1] * len(upper_triangular)
    for row in range(len(upper_triangular) - 1, -1, -1):

        # substitute in the known roots for all the non-zero coefficients and subtract from the right hand side for the current row
        for root in range(row + 1, len(upper_triangular)):
            c_vector[row] = c_vector[row] - roots[root] * upper_triangular[row][root]

        # calculate the root for this row with the updated right hand side value
        roots[row] = c_vector[row] / upper_triangular[row][row]

    return roots


def main():
    """ This will print the roots to the system of equations described by the left_hand_side matrix and the right hand side vector """

    left_hand_side = np.array([[4, 0, -2], [0, 1, 1], [-2, 1, 3]])
    right_hand_side = np.array([4, 2, 0])

    R, R_transpose = cholesky_factorization(left_hand_side)
    print solve_system_of_equations(R, R_transpose, right_hand_side)


if __name__ == "__main__":
    main()
