import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib import colors


def timestep(forest):
    temp = np.where(forest == 1, np.random.random(shape), tree)
    forest = np.where((0 <= temp) & (temp < ignition_probability), fire, forest)
    burning = np.where(forest == 2)
    print(burning)
    x = [burning[0][0], burning[1][0]]
    for x in range(0, len(burning[0])):
        coord = [burning[0][0], burning[1][0]]
        for diff in directions:
            if forest[coord[0]+diff[0], coord[1]+diff[1]] == tree:
                forest[coord[0] + diff[0], coord[1] + diff[1]] = fire
    return forest


directions = [(-1, -1), (-1, 0), (1, 0), (1, 1)]
ignition_probability = 0.01
growth_probability = 0.1
shape = (100, 100)
colour_list = colors.ListedColormap(['Black', 'Green', 'Red'])
grid = np.zeros(shape)
empty, tree, fire = 0, 1, 2

tempsetup = np.where(grid == empty, np.random.random(shape), empty)
grid = np.where((0 <= tempsetup) & (tempsetup < growth_probability), tree, grid)
grid = timestep(grid)


plt.imshow(grid, cmap=colour_list)
plt.show()
print(grid)


#plt.figure()

