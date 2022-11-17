import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib import colors

ignition_probability = 0.01
growth_probability = 0.1
shape = (100, 100)
colour_list = colors.ListedColormap(['Black', 'Green', 'Red'])
grid = np.zeros(shape)
empty, tree, fire = 0, 1, 2

temp = np.where(grid == empty, np.random.random(shape), empty)
grid = np.where((0 <= temp) & (temp < growth_probability), tree, grid)
temp = np.where(grid == 1, np.random.random(shape), tree)
grid = np.where((0 <= temp) & (temp < ignition_probability), fire, grid)


plt.imshow(grid, cmap=colour_list)
plt.show()
print(grid)


#plt.figure()

