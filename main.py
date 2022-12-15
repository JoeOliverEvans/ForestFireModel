import datetime
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib import colors
import matplotlib
import scipy.optimize as optimize
import scipy.stats as stats
import time
import os

# check if figures folder exists
if not os.path.isdir("figures"):
    # if the figures folder is not found create one
    os.makedirs("figures")


# TODO vary f and p values
# TODO calculate chisquared for fits


def timestep(forest, shape):
    """
    Advances the forest forward one iteration
    :param shape: shape of grid
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


def visual(shape):
    matplotlib.use("TkAgg")  # Necessary for my IDE to display animation
    colour_list = colors.ListedColormap(['Black', 'Green', 'Red'])
    grid = np.zeros(shape)  # Makes an empty forest
    grid[0, 0] = fire  # Necessary for the colormap to assign values to 0, 1 and 2 as they are based on the initial
    # maximum value doesn't impact the forest significantly

    fig, ax = plt.subplots()
    image = ax.imshow(grid, cmap=colour_list)

    def animate(frame):  # Generates the next image based on the forest
        image.set_data(animate.grid)
        animate.grid = timestep(animate.grid, shape)

    animate.grid = grid
    interval = 40  # Time between frames in ms
    amin = matplotlib.animation.FuncAnimation(fig, animate, interval=interval, frames=200)
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
                    labeled = np.where(labeled == left_label, above_label, labeled)  # Simple but inefficient merging
                    # of labels, used as union find failed to work as intended, not a problem for this investigation
    return labeled


def testing_Kopelman():
    """
    Algorithm to test Hoshen Kopelman implementation, using a 10x10 grid that contains "U" shapes, harder to label
    Alternatively it can also test a random 50% 1 50% 0 grid
    :return: None
    """
    grid = np.random.randint(0, 2, size=(10, 10), dtype=int)
    shape = (10, 10)
    hosh = Hoshen_Kopelman(grid, shape)
    plt.subplot(1, 2, 1)
    plt.title('Input Forest')
    plt.axis('off')
    image1 = plt.imshow(grid, cmap='seismic')
    for (j, i), label in np.ndenumerate(grid):
        plt.text(i, j, label, ha='center', va='center', color='white')
    plt.subplot(1, 2, 2)
    plt.title('Labelled Clusters')
    plt.axis('off')
    image2 = plt.imshow(hosh, cmap='seismic')
    for (j, i), label in np.ndenumerate(hosh):
        plt.text(i, j, label, ha='center', va='center', color='white')
    plt.tight_layout()
    plt.savefig('figures/Kopelman_test', dpi=240)
    plt.show()


def perimeter(labelled_grid, target, shape):
    """
    Returns the size of the cluster with number target
    :param shape: shape of forest
    :param labelled_grid: Grid resulting from Kopelman
    :param target: Number of cluster to find circumference of
    :return:
    """
    coords = np.argwhere(labelled_grid == target)
    counter = 0
    for coord in coords:
        minicounter = 0
        for direction in directions:
            if coord[0] + direction[0] < 0 or coord[0] + direction[0] >= shape[0] or coord[1] + direction[1] < 0 or \
                    coord[1] + direction[1] >= shape[1]:
                minicounter += 1
            else:
                if labelled_grid[coord[0] + direction[0], coord[1] + direction[1]] != target:
                    minicounter += 1
        counter += minicounter
    return counter


def perimeter_test():
    perimeter_test_grid = np.array([[0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                    [0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
                                    [0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
                                    [0, 0, 0, 0, 0, 1, 1, 0, 0, 0],
                                    [0, 0, 0, 0, 1, 1, 0, 0, 0, 0],
                                    [0, 0, 0, 1, 1, 0, 0, 0, 0, 0],
                                    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]])
    plt.subplot(1, 2, 1)  # Output in subplot to match font size etc. of other plots in report
    plt.title(f'Cluster perimeter : {perimeter(Hoshen_Kopelman(perimeter_test_grid, (10, 10)), 1, (10, 10))}')
    plt.axis('off')
    image1 = plt.imshow(perimeter_test_grid, cmap='seismic')
    for (j, i), label in np.ndenumerate(perimeter_test_grid):
        plt.text(i, j, label, ha='center', va='center', color='white')
    plt.tight_layout()
    plt.savefig('figures/perimeter_test', dpi=240)
    plt.show()


def powerlaw(x, k, a):
    return k * x ** -a


def powerfixed(x, k):
    return k * x ** (-2)


def investigation(shape, iterations, discard, animate, test):
    if test:
        testing_Kopelman()
    grid = np.zeros(shape)
    empties = [np.count_nonzero(grid == empty)]
    trees = [np.count_nonzero(grid == tree)]
    fires = [np.count_nonzero(grid == fire)]
    clusters, perimeters, radii = [], [], []
    t_0 = time.time()
    first = True
    sample = True
    s_0 = 0
    print(f'Time estimate will be given after discarded ({discard}) iterations')
    for x in range(1, iterations):
        grid = timestep(grid, shape)
        empties.append(np.count_nonzero(grid == empty))
        trees.append((np.count_nonzero(grid == tree)))
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
                    perimeters.append(perimeter(labelled_grid, y, shape))
                    clusters.append(int(np.count_nonzero(labelled_grid == y)))
                    # radii.append(calculateradius(labelled_grid, y, shape))
    xvalues = np.arange(0, iterations)
    rhobar = np.mean(trees[discard:]) / (shape[0]) ** 2
    plt.title(f'Number of trees in the forest over {iterations} iterations ' + r'$\rho ̄$' + f'$={rhobar:.3f}$')
    plt.errorbar(xvalues, empties, label='empties', fmt='.', markersize=4, color='k')
    plt.errorbar(xvalues, trees, label='trees', fmt='.', markersize=4, color='green')
    plt.xlabel('Number of iterations')
    plt.ylabel('Number of trees and empties')
    plt.legend()
    plt.savefig(f'figures/{shape}treenums.png', dpi=240)
    plt.show()

    datacolor = 'midnightblue'
    fitcolor = 'firebrick'
    bins = np.arange(1, np.max(clusters) + 1, 1, dtype=np.float64)
    cluster_hist_data = plt.hist(clusters, bins=bins, label='Cluster size', color=datacolor)
    cluster_yvals = cluster_hist_data[0]
    popt, pcov = optimize.curve_fit(powerlaw, xdata=cluster_hist_data[1][1:-1], ydata=cluster_yvals[1:],
                                    absolute_sigma=True)    # omitting cluster size 1
    topt, tcov = optimize.curve_fit(powerfixed, xdata=cluster_hist_data[1][1:-1], ydata=cluster_yvals[1:],
                                    absolute_sigma=True)    # fixed value of a=2
    cutoff5 = None
    for x in cluster_hist_data[1][:-1]:
        if powerlaw(x, *popt) < 5:
            cutoff5 = int(x - 1)
            break
    if cutoff5 is None:
        cutoff5 = int(len(bins) - 2)

    manualchi = np.sum((cluster_yvals[:cutoff5] - powerfixed(cluster_hist_data[1][:cutoff5], *topt)) ** 2 / powerfixed(
        cluster_hist_data[1][:cutoff5], *topt))
    manualchi1 = np.sum((cluster_yvals[:cutoff5] - powerlaw(cluster_hist_data[1][:cutoff5], *popt)) ** 2 / powerlaw(
        cluster_hist_data[1][:cutoff5], *popt))
    print(manualchi)
    print(manualchi1)
    print(stats.distributions.chi2.cdf(manualchi, cutoff5 - 1))
    # print(scipy.stats.chisquare(yvals[:cutoff5], powerlaw(y[1][:cutoff5], *popt)))

    print(topt)
    plt.plot(cluster_hist_data[1][:-1], powerlaw(cluster_hist_data[1][:-1], *popt),
             label=r'$kx^{-a}$' + f': $k = ${popt[0]:.2f}, $a = ${popt[1]:.2f}', color=fitcolor)
    plt.plot(cluster_hist_data[1][:-1], powerfixed(cluster_hist_data[1][:-1], *topt),
             label=r'$kx^{-2}$' + f': $k = ${topt[0]:.2f}', color='lightpink')
    plt.xlim(1, 30)
    plt.title(f'Frequency of cluster size for {shape} forest, {iterations} iterations')
    plt.ylabel('Frequency')
    plt.xlabel('Cluster size')
    plt.legend()
    plt.savefig(f'figures/{shape}area.png', dpi=240)
    plt.show()

    plt.title(f'Frequency of cluster size for {shape} forest, {iterations} iterations')
    plt.errorbar(cluster_hist_data[1][:-1], cluster_yvals, fmt='.', label='Cluster Size Frequency', color=datacolor)
    plt.plot(bins, powerlaw(bins, *popt), label=r'$kx^{-a}$' + f': $k = ${popt[0]:.2f}, $a = ${popt[1]:.2f}',
             color=fitcolor)
    plt.plot(cluster_hist_data[1][:-1], powerfixed(cluster_hist_data[1][:-1], *topt),
             label=r'$kx^{-2}$' + f': $k = ${topt[0]:.2f}', color='lightpink')
    plt.ylabel('Frequency')
    plt.xlabel('Cluster size')
    plt.loglog()
    plt.legend()
    plt.savefig(f'figures/{shape}logarea.png', dpi=240)
    plt.show()

    bincirc = np.arange(4, np.max(perimeters) + 2, 2, dtype=np.float64)
    perimeter_hist_data = plt.hist(perimeters, bins=bincirc, label='Circumference', color=datacolor)
    perimeter_yvals = perimeter_hist_data[0]
    popt, pcov = optimize.curve_fit(powerlaw, xdata=perimeter_hist_data[1][1:-1], ydata=perimeter_yvals[1:],
                                    absolute_sigma=True)

    plt.plot(perimeter_hist_data[1][:-1], powerlaw(perimeter_hist_data[1][:-1], *popt),
             label=r'$kx^{-a}$' + f': $k = ${popt[0]:.2f}, $a = ${popt[1]:.2f}', color=fitcolor)
    plt.xlim(1, 30)
    plt.title(f'Frequency of perimeter fitted with a power law for {shape}')
    plt.ylabel('Frequency')
    plt.xlabel('Perimeter')
    plt.legend()
    plt.savefig(f'figures/{shape}circ.png', dpi=240)
    plt.show()

    plt.title(f'Frequency of perimeter fitted with a power law for {shape}')
    plt.errorbar(perimeter_hist_data[1][:-1], perimeter_yvals, fmt='.', label='Cluster Size Frequency', color=datacolor)
    plt.plot(bincirc, powerlaw(bincirc, *popt), label=r'$kx^{-a}$' + f': $k = ${popt[0]:.2f}, $a = ${popt[1]:.2f}',
             color=fitcolor)
    plt.ylabel('Frequency')
    plt.xlabel('Perimeter')
    plt.loglog()
    plt.legend()
    plt.savefig(f'figures/{shape}logcirc.png', dpi=240)
    plt.show()
    return


def Gigafunction(shapes=None, iterations=3000, discard=100, animate=False, test=False):
    if shapes is None:
        shapes = [(50, 50), (200, 200)]
    for shape in shapes:
        investigation(shape, iterations, discard, animate, test)
    if animate:
        visual(shape=(100, 100))


directions = [(0, -1), (-1, 0), (1, 0), (0, 1)]  # Used to index neighboring cells in the forest
growth_probability = 0.01
ignition_probability = growth_probability * 1 / 70
empty, tree, fire = 0, 1, 2  # Assigning values to represent empty, tree and fire
time_sample_iterations = 10  # Used for the estimated completion time

# visual()
# study_numbers()
# testing_Kopelman()
# investigating_clusters(200)
# investigation((50, 50), 500, 100, False, False)
Gigafunction(shapes=[(50, 50)], iterations=10000, animate=False, discard=9000)
# visual((100, 100))
# perimeter_test()
