import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib import colors
import matplotlib


def timestep(forest):
    """
    Advances the forest forward one iteration
    :param forest: grid containing trees and fire
    :return: new grid propagated due to the rules
    """
    fires = np.argwhere(forest == fire)
    random = np.random.random(shape)
    new_forest = np.where(forest == tree, tree, 0)
    new_growth = np.where((forest == empty) & (random < growth_probability), tree, 0)
    ignition = np.where((forest == tree) & (random < ignition_probability), 1, 0)
    new_forest = new_forest + new_growth + ignition
    for x in fires:
        for y in directions:
            try:
                if forest[x[0] + y[0]][x[1] + y[1]] == tree:
                    new_forest[x[0] + y[0]][x[1] + y[1]] = fire
            except IndexError:
                pass
    return new_forest


directions = [(0, -1), (-1, 0), (1, 0), (0, 1)]
ignition_probability = 0.001
growth_probability = 0.01
shape = (100, 100)
empty, tree, fire = 0, 1, 2


def visual():
    matplotlib.use("TkAgg")
    colour_list = colors.ListedColormap(['Black', 'Green', 'Red'])
    grid = np.zeros(shape)
    for p in range(0, 50):
        grid = timestep(grid)

    fig, ax = plt.subplots()
    image = ax.imshow(grid, cmap=colour_list)

    def animate(frame):
        image.set_data(animate.grid)
        animate.grid = timestep(animate.grid)

    animate.grid = grid
    interval = 50
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
    ax1.errorbar(xvalues, fires, label='fires', fmt='.', color='red')
    ax1.set_xlabel('number of iterations')

    plt.legend()
    fig.tight_layout()
    plt.show()


def Hoshen_Kopelman(grid):
    """
    Python version of algorithm written in pseudocode on https://en.wikipedia.org/wiki/Hoshen%E2%80%93Kopelman_algorithm
    Had to make significant changes to check for clusters that were adjacent with different labels
    :param grid: Forest containing empty, trees and fires
    :return:numbered grid of clusters
    """
    labeled = np.zeros(shape, dtype=int)
    largest_label = int(0)
    equivalent = []
    for xkop in range(0, shape[0]):
        for ykop in range(0, shape[1]):
            if grid[xkop, ykop] == 1:
                if ykop - 1 < 0:
                    left = 0
                else:
                    left = grid[xkop, ykop - 1]
                if xkop - 1 < 0:
                    above = 0
                else:
                    above = grid[xkop - 1, ykop]
                if left != tree and above != tree:
                    largest_label += 1
                    labeled[xkop, ykop] = largest_label
                if left == tree and above != tree:
                    labeled[xkop, ykop] = labeled[xkop, ykop - 1]
                if left != tree and above == tree:
                    labeled[xkop, ykop] = labeled[xkop - 1, ykop]
                if left == tree and above == tree:
                    labeled[xkop, ykop] = labeled[xkop, ykop - 1]
                    if labeled[xkop - 1, ykop] != labeled[xkop, ykop - 1]:
                        join = [labeled[xkop - 1, ykop], labeled[xkop, ykop - 1]]
                        join.sort()
                        if join not in equivalent:
                            equivalent.append(join)

    done = False
    while not done:
        change = False
        for n in range(0, len(equivalent)):
            for m in range(0, len(equivalent)):
                if n != m:
                    for x in equivalent[n]:
                        if x in equivalent[m]:
                            equivalent[n] = [*{*equivalent[n], *equivalent[m]}]
                            equivalent.pop(m)
                            change = True
                        if change:
                            break
                if change:
                    break
            if change:
                break
        if not change:
            done = True

    for pair in equivalent:
        pair.sort()
        for r in range(1, len(pair)):
            labeled = np.where(labeled == pair[r], pair[0], labeled)
    p = 1
    while p < largest_label:
        if p not in labeled and p < np.max(labeled):
            labeled = np.where(labeled > p, labeled - 1, labeled)
        else:
            p += 1
    return labeled


def testing_Hoshen():
    """
    Algorithm to test Hoshen Kopelman implementation, success and fail messages are regarding all labels being in
    order not the accuracy of labelling the clusters.
    :return: None
    """
    grid = np.zeros(shape)

    for x in range(0, 50):
        grid = timestep(grid)
    print(grid)
    hosh = Hoshen_Kopelman(grid)
    print(hosh)
    new = np.where(hosh != 0, 1, 0)
    ma = np.max(hosh)
    fail = False
    for x in range(0, ma + 1):
        if x not in hosh:
            fail = True

    if not fail:
        print('success')

    else:
        print('fail')
    return None


def avg_cluster_size(labelled_grid):
    """
    Returns the average cluster size of groups of trees
    :param labelled_grid: forest with Kopelman labels
    :return: avg: average size of tree clusters
    """
    runningtot = 0
    for x in range(1, np.max(labelled_grid) + 1):
        runningtot += np.count_nonzero(labelled_grid == x)
    avg = runningtot / np.max(labelled_grid)
    return avg


def investigating_clusters(iterations=5*10**2, repeats=20):
    """
    Records number of clusters
    :param iterations:
    :param repeats:
    :return:
    """
    number_of_clusters = []
    average_size_of_cluster = []
    grid = np.zeros(shape)
    for y in range(0, repeats):
        for x in range(0, iterations):
            grid = timestep(grid)
        labelled_grid = Hoshen_Kopelman(grid)
        number_of_clusters.append(np.max(labelled_grid))
        average_size_of_cluster.append(avg_cluster_size(labelled_grid))
    xvalues = np.linspace(iterations, iterations*repeats, repeats, dtype=int)
    plt.errorbar(xvalues, number_of_clusters)
    plt.show()
    plt.errorbar(xvalues, average_size_of_cluster)
    plt.show()
    print(number_of_clusters)
    print(average_size_of_cluster)


#study_numbers()
#testing_Hoshen()
investigating_clusters()

