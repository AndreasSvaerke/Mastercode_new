import numpy as np
import matplotlib.pyplot as plt

def plot(filename, epsilon):
    mean = np.nanmean(filename)
    std = np.nanstd(filename)

    plt.imshow(filename, vmin = mean - epsilon * std , vmax = mean + epsilon * std)


def plot_transient_phy(data, positions, eps = 1, inte = 25):
    x = positions[0]
    y = positions[1]
    plot(data[int(y-inte):int(y+inte) , int(x-inte):int(x+inte)] , eps)
    plt.show()
