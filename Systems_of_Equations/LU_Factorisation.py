"""
Script uses PA=LU factorisation to solve a system of equations
This is the established method of using Guassina elimination to solve a system of equations
The method works by converting the system of equations in matirx form from Ax = b to PLUx = Pb

Here P is a pivot matrix which represents the row exchanges made during row redution with partial pivoting
L is a lower trangular matrix which represents the multipiers used in the row redutcion
U is an upper triangular matrix which is the result of the row reductions

we can then solve the system of equations by first solving (1) for c and the (2) for x where c is a vector
(1) Lc = Pb
(2) Ux = c
"""

import numpy as np


def swap_rows(matrix, row_1, row_2):
    """
    Swaps 2 rows in a matrix

    matrix: (numpy 2D array), Matrix containing rows to be swapped
    row_1: (int), Index of the first row to be swapped, indexing starts at 0
    row_1: (int), Index of the second row to be swapped, indexing starts at 0
    """

    matrix[[row_1, row_2]] = matrix[[row_2, row_1]]


def lu_factorisation(left_hand_side):
    """ Performs PA = LU factorisation, This is the complete form of Guassian elimination

    Row exchanges are used to solve the issues of zeros at partial pivots and swamping
    These row exchanges involve choosing the largest possible value for the pivot points
    The matrix P is used to keep track of the row exchanges

    left_hand_side: (numpy 2D array), matrix represnting a system of equations to be factorised
    """

    p_matrix = np.identity(len(left_hand_side))
    l_matrix = np.zeros((len(left_hand_side), len(left_hand_side)))
    u_matrix = np.copy(left_hand_side)

    # perform row operations and partial pivoting to get upper and lower triangular matricies
    for row_top in range(len(u_matrix)):

        # find row with largest value in pivot column and swap with pivot row
        max_column_value = np.abs(u_matrix[row_top][row_top])
        max_column_index = row_top

        for row_bottom in range(row_top+1, len(u_matrix)):
            if np.abs(u_matrix[row_bottom][row_top]) > max_column_value:
                max_column_value = np.abs(u_matrix[row_bottom][row_top])
                max_column_index = row_bottom

        # swap the pivot rows
        swap_rows(u_matrix, row_top, max_column_index)
        swap_rows(l_matrix, row_top, max_column_index)
        swap_rows(p_matrix, row_top, max_column_index)

        # perform row operations to reduce matrix to upper triangular
        for row_bottom in range(row_top+1, len(u_matrix)):

            # check for zero pivots
            if (u_matrix[row_top][row_top]) == 0:
                raise ValueError("There is a zero pivot point at %s, %s" % (row_top, row_top))

            # calculate the value to multiply the top row by to remove the bottom row
            multiplier = float(u_matrix[row_bottom][row_top]) / float(u_matrix[row_top][row_top])

            # subtract the top row form the bottom row (left hand side)
            for element in range(row_top, len(u_matrix)):
                u_matrix[row_bottom][element] = u_matrix[row_bottom][element] - (multiplier * u_matrix[row_top][element])

            # add the multiplier into the corresponding position in the lower triangular matrix
            l_matrix[row_bottom][row_top] = multiplier

    l_matrix += np.identity(len(left_hand_side))

    return p_matrix, l_matrix, u_matrix


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

    # declare the system of equations as a matrix (left hand side) and a vector (right hand side)
    left_hand_side = np.array([[2, 1, 5], [4, 4, -4], [1, 3, 1]])
    right_hand_side = np.array([5, 0, 6])

    p_matrix, l_matrix, u_matrix = lu_factorisation(left_hand_side)

    print solve_system_of_equations(u_matrix, l_matrix, p_matrix.dot(right_hand_side))

if __name__ == "__main__":
    main()
