"""
Interpolation is the reverse of evalutaion, instead of being given a polynomial and finding points on it
We are given points and asked to find the polynomial that the points exist on

We are guaranteed that a given N points there is a polynomial order at most N - 1 that interplates these points
The proof of this is that we need N - 1 turning points to be certain to reach all data points
Increasing the order of the polynomial adds an extra tunring point

Lagrange Interpolation indroduces a formula which creates a N - 1 polinomial for N points
i.e. If we had three points (x1, y1), (x2, y2) and (x3, y3)

P(x) =      y1 * (x - x2)*(x - x3) / (x1 - x2)*(x1 - x3)
        +   y2 * (x - x1)*(x - x3) / (x2 - x1)*(x2 - x3)
        +   y3 * (x - x1)*(x - x2) / (x3 - x1)*(x3 - x2)

From this formula we can see that when x is substituted with x1 the polinomial evaluates to y1
The same is ture for the (x2, y2) and (x3, y3) pairs
"""


import numpy as np
from matplotlib import pyplot as plt


def plot_interpolations(interpolation_points, polynomial_points):
    """ plots the computed polynomial and the interpolation points used to calculate it """

    x_poly, y_poly = zip(*polynomial_points)
    x_interp, y_interp = zip(*interpolation_points)

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)

    # center the axis at (0, 0)
    ax.spines['left'].set_position('center')
    ax.spines['bottom'].set_position('center')

    plt.plot(x_poly, y_poly, color='blue', linewidth=2)                      # plot polynomial
    plt.plot(x_interp, y_interp, 'bo', color='red', linewidth=2)             # plot interpolation points

    # add lables, legend and title
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend(['f1(x)', "f2(x)"], loc=2)
    plt.title("Multivariate Newton Raphson Method")

    # show the plot
    plt.show()


def lagrange_interpolation(interpolation_points, test_points):
    """
    Uses lagrange interpolation to determine a polinomial solution
    given a set of 2D points which must lie on the polinomial

    interpolation_points: (numpy 2D array) List of points
    test_points: (numpy 1D array) List of x-values the polinomial is to be evaluated at
    """

    def factor(index):
        """ Calculates the fractions that's mulitpled by each y-value """

        f = [(float(x) - interpolation_points[i][0]) / (interpolation_points[index][0] - interpolation_points[i][0])
             for i in range(num_points) if i != index]

        # multiply to list elements to from the complete fraction
        return np.prod(np.array(f))


    results = list()                                # will store the x and y values for the test_points evaluated on the polynomial
    num_points = len(interpolation_points)          # Number of interpolation points

    # for each test point calculate the y value
    for x in test_points:
        y = np.sum([factor(j) * interpolation_points[j][1] for j in range(num_points)])
        results.append(np.array([x, y]))

    return np.array(results)

def main():

    interpolation_points = np.array([[0, 1], [2, 2], [3, 4]])
    test_points = np.arange(-10.0, 10.0, 0.1)

    results = lagrange_interpolation(interpolation_points, test_points)
    plot_interpolations(interpolation_points, results)

if __name__ == "__main__":
    main()
