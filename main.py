import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib import colors
import matplotlib

matplotlib.use("TkAgg")


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
    ignition = np.where((forest == tree) & (random < ignition_probability), tree, 0)
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
colour_list = colors.ListedColormap(['Black', 'Green', 'Red'])
grid = np.zeros(shape)
empty, tree, fire = 0, 1, 2


for p in range(0, 100):
    grid = timestep(grid)

fig, ax = plt.subplots()
image = ax.imshow(grid, cmap=colour_list)


def animate(frame):
    image.set_data(animate.grid)
    animate.grid = timestep(animate.grid)


animate.grid = grid
interval = 50
animation = matplotlib.animation.FuncAnimation(fig, animate, interval=interval, frames=200)

plt.show()
