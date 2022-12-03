import datetime

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib import colors
import matplotlib
import scipy.optimize
import time


def timestep(forest):
    """
    Advances the forest forward one iteration
    :param forest: grid containing trees and fire
    :return: new grid propagated due to the rules
    """
    fires = np.argwhere(forest == fire)  # Find the location of fires
    random = np.random.random(shape)  # Generate random values for the forest to be used for new growth etc.
    new_forest = np.where(forest == tree, tree, 0)
    new_growth = np.where((forest == empty) & (random < growth_probability), tree, 0)  # Grows new trees
    ignition = np.where((forest == tree) & (random < ignition_probability), 1, 0)  # Fire value is 2 however we will
    # add these to trees, so it will become 2 when added
    new_forest = new_forest + new_growth + ignition  # Constructs the new forest
    for x in fires:  # Propagates the fire from the previous forest
        for y in directions:
            try:  # Catches when the fire spreads out of bounds
                if forest[x[0] + y[0]][x[1] + y[1]] == tree:
                    new_forest[x[0] + y[0]][x[1] + y[1]] = fire
            except IndexError:
                pass
    return new_forest


def visual():
    matplotlib.use("TkAgg")  # Necessary for my IDE to display animation
    colour_list = colors.ListedColormap(['Black', 'Green', 'Red'])
    grid = np.zeros(shape)  # Makes an empty forest
    grid[0, 0] = fire  # Necessary for the colormap to assign values to 0, 1 and 2 as they are based on the initial
    # maximum value doesn't impact the forest significantly

    fig, ax = plt.subplots()
    image = ax.imshow(grid, cmap=colour_list)

    def animate(frame):  # Generates the next image based on the forest
        image.set_data(animate.grid)
        animate.grid = timestep(animate.grid)

    animate.grid = grid
    interval = 50  # Time between frames in ms
    amin = matplotlib.animation.FuncAnimation(fig, animate, interval=interval, frames=200)
    plt.show()


def study_numbers(n=10 ** 3):
    """
    Analyse the number of trees fires and empty squares over time
    :param: number of iterations
    :return: Nothing
    """
    grid = np.zeros(shape)
    empties = [np.count_nonzero(grid == empty)]
    trees = [np.count_nonzero(grid == tree)]
    fires = [np.count_nonzero(grid == fire)]
    for x in range(1, n):
        grid = timestep(grid)
        empties.append(np.count_nonzero(grid == empty))
        trees.append((np.count_nonzero(grid == tree)))
        fires.append((np.count_nonzero(grid == fire)))

    xvalues = np.arange(0, n)

    fig, ax1 = plt.subplots()

    ax1.errorbar(xvalues, empties, label='empties', fmt='.', color='k')
    ax1.errorbar(xvalues, trees, label='trees', fmt='.', color='green')
    ax1.set_xlabel('number of iterations')
    ax1.set_ylabel('number of trees and empties')

    plt.legend()
    fig.tight_layout()
    plt.show()


def Hoshen_Kopelman(grid, shape):
    """
    Python version of algorithm written in pseudocode on https://en.wikipedia.org/wiki/Hoshen%E2%80%93Kopelman_algorithm
    Had to make significant changes to check for clusters that were adjacent with different labels
    :param grid: Forest containing empty, trees and fires
    :return:numbered grid of clusters
    """
    labeled = np.zeros(shape, dtype=int)  # Creates the forest
    largest_label = int(0)  # Variable to track largest label
    for xkop in range(0, shape[0]):  # Iterate through the grid
        for ykop in range(0, shape[1]):
            if grid[xkop, ykop] == tree:
                if ykop - 1 < 0:  # Check if left or above is in or out of the array
                    left = 0
                else:
                    left = grid[xkop, ykop - 1]
                if xkop - 1 < 0:
                    above = 0
                else:
                    above = grid[xkop - 1, ykop]
                if left != tree and above != tree:  # Start a new label as there are no bordering trees
                    largest_label += 1
                    labeled[xkop, ykop] = largest_label
                elif left == tree and above != tree:  # Assign label from the left
                    labeled[xkop, ykop] = labeled[xkop, ykop - 1]
                elif left != tree and above == tree:  # Assign label from the right
                    labeled[xkop, ykop] = labeled[xkop - 1, ykop]
                else:  # Assign label from the top and change all labels to the left to above label
                    labeled[xkop, ykop] = labeled[xkop - 1, ykop]
                    left_label = labeled[xkop, ykop - 1]
                    above_label = labeled[xkop - 1, ykop]
                    labeled = np.where(labeled == left_label, above_label, labeled)     # Simple but inefficient
                    # merging of labels
    return labeled


def testing_Kopelman():
    """
    Algorithm to test Hoshen Kopelman implementation, using a 10x10 grid that contains "U" shapes, harder to label
    :return: None
    """
    # TODO make this look good
    shape = (10, 10)
    grid = np.array([[0, 1, 0, 0, 0, 0, 1, 1, 1, 0],
                     [0, 0, 0, 0, 1, 1, 0, 1, 0, 1],
                     [0, 0, 1, 0, 0, 1, 1, 1, 1, 1],
                     [0, 1, 1, 1, 1, 0, 1, 0, 0, 0],
                     [0, 0, 0, 0, 1, 1, 1, 0, 0, 0],
                     [0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
                     [0, 1, 0, 1, 0, 0, 1, 0, 1, 1],
                     [0, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                     [0, 0, 0, 0, 0, 1, 0, 1, 1, 1],
                     [0, 0, 0, 0, 0, 0, 0, 0, 1, 0]])
    print(grid)
    hosh = Hoshen_Kopelman(grid, shape)
    print(hosh)


def powerlaw(x, k, a):
    return k * x ** a


def investigating_clusters(iterations=3 * 10 ** 3, discard=100):
    """
    Records frequency of different sized clusters
    :param discard:
    :param iterations:
    :return:
    """
    clusters = []
    grid = np.zeros(shape)
    t_0 = time.time()
    first = True
    sample = True
    s_0 = 0
    print(f'Time estimate will be given after discarded ({discard}) iterations')
    for x in range(0, iterations):
        grid = timestep(grid)
        if (time.time() - t_0) > 5:
            print(str(x) + ' / ' + str(iterations))
            t_0 = time.time()
        if x >= discard:
            if sample:
                sample = False
                s_0 = time.time()
            if x >= (discard + time_sample_iterations):
                if first:
                    time_to_completion = round(((time.time() - s_0) * (iterations - discard) / time_sample_iterations))
                    print('Estimated time to completion : ' + str(datetime.timedelta(seconds=time_to_completion)))
                    first = False
            labelled_grid = Hoshen_Kopelman(grid, shape)
            for y in range(1, np.max(labelled_grid) + 1):
                if y in labelled_grid:  # Avoid appending 0 values for missing labels
                    clusters.append(int(np.count_nonzero(labelled_grid == y)))
    bins = np.arange(1, np.max(clusters) + 1, 1)
    y = plt.hist(clusters, bins=bins, label='Cluster size')
    yvals = y[0]
    popt, pcov = scipy.optimize.curve_fit(powerlaw, xdata=bins[:-1], ydata=yvals)
    plt.plot(bins, powerlaw(bins, *popt), label=fr'$kx^a$ : $k = ${popt[0]:.2f}, $a = ${popt[1]:.2f}')
    plt.xlim(1, 30)
    plt.title(f'Frequency of cluster size fitted with a power law for {shape}')
    plt.ylabel('Frequency')
    plt.xlabel('Clusters')
    plt.legend()
    plt.show()


def investigation(shape, iterations, discard, animate, test):
    if test:
        testing_Kopelman()
    grid = np.zeros(shape)
    empties = [np.count_nonzero(grid == empty)]
    trees = [np.count_nonzero(grid == tree)]
    fires = [np.count_nonzero(grid == fire)]
    clusters = []
    t_0 = time.time()
    first = True
    sample = True
    s_0 = 0
    print(f'Time estimate will be given after discarded ({discard}) iterations')
    for x in range(0, iterations):
        grid = timestep(grid)
        empties.append(np.count_nonzero(grid == empty))
        trees.append((np.count_nonzero(grid == tree)))
        fires.append((np.count_nonzero(grid == fire)))
        if (time.time() - t_0) > 5:
            print(str(x) + ' / ' + str(iterations))
            t_0 = time.time()
        if x >= discard:
            if sample:
                sample = False
                s_0 = time.time()
            if x >= (discard + time_sample_iterations):
                if first:
                    time_to_completion = round(((time.time() - s_0) * (iterations - discard) / time_sample_iterations))
                    print('Estimated time to completion : ' + str(datetime.timedelta(seconds=time_to_completion)))
                    first = False
            labelled_grid = Hoshen_Kopelman(grid)
            for y in range(1, np.max(labelled_grid) + 1):
                if y in labelled_grid:  # Avoid appending 0 values for missing labels
                    clusters.append(int(np.count_nonzero(labelled_grid == y)))

    xvalues = np.arange(0, iterations)
    plt.errorbar(xvalues, empties, label='empties', fmt='.', color='k')
    plt.errorbar(xvalues, trees, label='trees', fmt='.', color='green')
    plt.xlabel('Number of iterations')
    plt.ylabel('Number of trees and empties')
    plt.legend()
    plt.show()

    bins = np.arange(1, np.max(clusters) + 1, 1)
    y = plt.hist(clusters, bins=bins, label='Cluster size')
    yvals = y[0]
    popt, pcov = scipy.optimize.curve_fit(powerlaw, xdata=bins[:-1], ydata=yvals)
    plt.plot(bins, powerlaw(bins, *popt), label=fr'$kx^a$ : $k = ${popt[0]:.2f}, $a = ${popt[1]:.2f}')
    plt.xlim(1, 30)
    plt.title(f'Frequency of cluster size fitted with a power law for {shape}')
    plt.ylabel('Frequency')
    plt.xlabel('Clusters')
    plt.legend()
    plt.show()
    return


def Gigafunction(shapes=None, iterations=2000, discard=100, animate=False, test=False):
    if shapes is None:
        shapes = [(50, 50), (200, 200)]
    if animate:
        visual()
    for shape in shapes:
        investigation(shape, iterations, discard, animate, test)


directions = [(0, -1), (-1, 0), (1, 0), (0, 1)]
ignition_probability = 0.001
growth_probability = 0.01
shape = (50, 50)
empty, tree, fire = 0, 1, 2
time_sample_iterations = 10

# visual()
# study_numbers()
testing_Kopelman()
# investigating_clusters(2000)
# Gigafunction()
