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
shape = (50, 50)
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
    label = np.zeros(shape, dtype=int)
    labels = np.arange(0, (shape[0]*shape[1])/2)
    for xkop in range(0, shape[0]):
        for ykop in range(0, shape[1]):
            if grid[xkop, ykop] == tree:
                if ykop-1<0:
                    left = 0
                else:
                    left = grid[xkop, ykop-1]
                if xkop-1<0:
                    above=0
                else:
                    above = grid[xkop-1, ykop]
                if left != tree and above != tree:
                    largest_label = largest_label + 1
                    label[xkop, ykop] = largest_label
                elif left == tree and above != tree:
                    label[xkop, ykop] = find(left)
                elif left != tree and above == tree:
                    label[xkop, ykop] = find(above)
                else:
                    union(left, above)
                    label[xkop, ykop] = find(left)
    '''for (int i=0; i < m; i++)
        for (int j=0; j < n; j++)
            if (matrix[i][j]) {
            int x = uf_find(matrix[i][j]);
            if (new_labels[x] == 0) {
            new_labels[0]++;
            new_labels[x] = new_labels[0];
            }
            matrix[i][j] = new_labels[x];
            }

    int
    total_clusters = new_labels[0];

    free(new_labels);'''
    new_labels = np.arange(0,max(labels))
    for xkop in range(0, shape[0]):
        for ykop in range(0, shape[1]):
            if grid[xkop, ykop] != 0:
                x = find(label[xkop,ykop])
                if new_labels[x] == 0:
                    new_labels[0] += 1
                    new_labels[x] = new_labels[0]
                label[xkop,ykop] = new_labels[x]
    return label


def Hoshen_Kopelman2(grid):
    labeled = np.zeros(shape, dtype=int)
    largest_label = int(0)
    equivalent = []
    for xkop in range(0, shape[0]):
        for ykop in range(0, shape[1]):
            if grid[xkop, ykop] == 1:
                if ykop-1<0:
                    left = 0
                else:
                    left = grid[xkop, ykop-1]
                if xkop-1<0:
                    above=0
                else:
                    above = grid[xkop-1, ykop]
                if left != tree and above != tree:
                    largest_label += 1
                    labeled[xkop, ykop] = largest_label
                if left == tree and above != tree:
                    labeled[xkop, ykop] = labeled[xkop, ykop-1]
                if left != tree and above == tree:
                    labeled[xkop, ykop] = labeled[xkop-1, ykop]
                if left == tree and above == tree:
                    labeled[xkop, ykop] = labeled[xkop, ykop-1]
                    if labeled[xkop-1, ykop] != labeled[xkop, ykop-1]:
                        join = [labeled[xkop-1, ykop], labeled[xkop, ykop-1]]
                        new = True
                        for n, h in enumerate(equivalent):
                            if join[0] in h and join[1] not in h:
                                equivalent[n].append(join[1])
                                new = False
                            elif join[1] in h and join[0] not in h:
                                equivalent[n].append(join[0])
                                new = False
                        if new:
                            equivalent.append(join)
    print(equivalent)
    for pair in equivalent:
        pair.sort()
        for r in range(1, len(pair)):
            labeled = np.where(labeled == pair[r], pair[0], labeled)
    print(labeled)
    for p in range(0, largest_label):
        if p not in labeled:
            np.where(labeled>p, labeled-1, labeled)
    '''labeled[xkop, ykop] = labeled[xkop-1, ykop]
    cond = True
    p = 0
    while cond:
        p += 1
        if labeled[xkop, ykop-1] != 0:
            labeled[xkop, ykop-p] = labeled[xkop, ykop]
        cond = False'''
    return labeled



#study_numbers()
grid = np.zeros(shape)
for x in range(0, 50):
    grid = timestep(grid)

print(grid)
hosh = Hoshen_Kopelman2(grid)
print(hosh)
new = np.where(hosh != 0, 1, 0)
print(np.max(np.where(grid==1, 1, 0)-new))



