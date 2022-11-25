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
shape = (10, 10)
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


def study_numbers(n=10**3):
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
    Python version of algoithm detailed in psuedocode on https://en.wikipedia.org/wiki/Hoshen%E2%80%93Kopelman_algorithm
    :param grid:
    :return:
    """
    def union(a, b):
        labels[find(a)] = find(b)

    def find(p):
        q = p
        while labels[q] != q:
            q = labels[q]
        while labels[p] != p:
            r = labels[p]
            labels[p] = q
            p = r
        return q

    largest_label = 0
    label = np.zeros(shape)
    labels = np.arange(0, shape[0]*shape[1])
    for xkop in range(0, shape[0]):
        for ykop in range(0, shape[1]):
            if grid[xkop, ykop] == tree:
                left = grid[xkop - 1, ykop]
                above = grid[xkop, ykop - 1]
                if left != tree and above != tree:
                    largest_label += 1
                    label[xkop, ykop] = largest_label
                elif left == tree and above != tree:
                    label[xkop, ykop] = find(left)
                elif left != tree and above == tree:
                    label[xkop, ykop] = find(above)
                else:
                    union(left, above)
                    label[xkop, ykop] = find(left)
    return label


def Hoshen_Kopelman2():


#study_numbers()
grid = np.zeros(shape)
for x in range(0, 50):
    grid = timestep(grid)

print(grid)
print(Hoshen_Kopelman2(grid))
