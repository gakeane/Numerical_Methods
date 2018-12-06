"""
This script uses guassian elimination to solve a system of n equations with n unknonwns
Guassian elimination consists of three operations

(1) Swap one equation for another.
(2) Add or subtract a multiple of one equation from another.
(3) Multiply an equation by a nonzero constant.

The goal is to use these operations to convert the matrix into upper triangular form (i.e. the bottom right triangle is all zeros)
We can then use back substitution (subbing in variables once we work out what they are to solve the system of equations)

The system of equations should consts of an NxN matrix representing the right hand side of the equaltions and a
Nx1 vector representing the left hand side of the equations

Note that the diagonal elements of the matrix are called pivots
If we have a zero pivot (i.e a_(ii) is zero) then the naive version of guassian elimination will halt (division by zero)
We need to swap the row with the zero pivot point woth the (largest) row below it -> partial pivoting
Partial pivoting invloves swaping the current row with the row that has the largest value in the current column
i.e. swap row i with row j if a_ij >= a_ik for all 1 <= k <= n

In the case where all entries in a column is zero the matrix is singular
Guassian elimination will always fail if the left hand side matrix is singular (non-invertable)
"""

import numpy as np


def guassian_elimination(left_hand_side, right_hand_side):
    """ Performs guassian elimination to calculate the roots of a system of equations

    left_hand_side: (2D numpy array matrix containing the left hand side of the system of equations
    right_hand_side: (1D numpy array) vector containing the right hand side of the vector

    NOTE: this version of the function hasn't been implemented to deal with zeros at pivot points
    To fix this when we have a zero at a pivot point we need swap the current row with one lower (in the matrix) than it
    """

    left_hand_shape = left_hand_side.shape
    right_hand_shape = right_hand_side.shape

    if len(left_hand_shape) != 2:
        raise ValueError("left hand side must be a 2D matrix, current dimensions: %s" % len(left_hand_shape))

    if len(right_hand_shape) != 1:
        raise ValueError("right hand side must be a 1D vector, current dimensions: %s" % len(right_hand_shape))

    if left_hand_shape[0] != left_hand_shape[1]:
        raise ValueError("left hand side must be a square matrix, current shape: %s" % left_hand_shape)

    if left_hand_shape[1] != len(right_hand_side):
        raise ValueError("left hand side and right hand side have incompatible dimesions: %s and %s" % (left_hand_shape, right_hand_shape))

    # perform row operations to reduce matrix to upper triangular
    for row_top in range(len(left_hand_side)):
        for row_bottom in range(row_top+1, len(left_hand_side)):

            # check for zero pivots
            if (left_hand_side[row_top][row_top]) == 0:
                raise ValueError("There is a zero pivot point at %s, %s" % (row_top, row_top))

            # calculate the value to multiply the top row by to remove the bottom row
            multiplier = float(left_hand_side[row_bottom][row_top]) / float(left_hand_side[row_top][row_top])

            # subtract the top row form the bottom row (left hand side)
            for element in range(row_top, len(left_hand_side)):
                left_hand_side[row_bottom][element] = left_hand_side[row_bottom][element] - (multiplier * left_hand_side[row_top][element])

            # subtract the top row from the bottom row (right hand side)
            right_hand_side[row_bottom] = right_hand_side[row_bottom] - multiplier * right_hand_side[row_top]


    # use back substitution to solve for the unknowns
    roots = [1] * len(left_hand_side)
    for row in range(len(left_hand_side) - 1, -1, -1):

        # substitute in the known roots for all the non-zero coefficients and subtract from the right hand side for the current row
        for root in range(row + 1, len(left_hand_side)):
            right_hand_side[row] = right_hand_side[row] - roots[root] * left_hand_side[row][root]

        # calculate the root for this row with the updated right hand side value
        roots[row] = right_hand_side[row] / left_hand_side[row][row]

    return roots


def main():

    # declare the system of equations as a matrix (left hand side) and a vector (right hand side)
    left_hand_side = np.array([[1, 2, -1], [2, 1, -2], [-3, 1, 1]])
    right_hand_side = np.array([3, 3, -6])

    left_hand_side = np.array([[-3, 1], [2, 4]])
    right_hand_side = np.array([-1, -4])

    print guassian_elimination(left_hand_side, right_hand_side)


if __name__ == "__main__":
    main()
