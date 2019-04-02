import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import os
import pandas as pd
from scipy.spatial import ConvexHull
import seaborn as sns
import nestle
########################
#Functions
##########################
def plot_ellipsoid_3d(ell, ax):
    """Plot the 3-d Ellipsoid ell on the Axes3D ax."""

    # points on unit sphere
    u = np.linspace(0.0, 2.0 * np.pi, 100)
    v = np.linspace(0.0, np.pi, 100)
    z = np.outer(np.cos(u), np.sin(v))
    y = np.outer(np.sin(u), np.sin(v))
    x = np.outer(np.ones_like(u), np.cos(v))

    # transform points to ellipsoid
    for i in range(len(x)):
        for j in range(len(x)):
            x[i,j], y[i,j], z[i,j] = ell.ctr + np.dot(ell.axes,
                                                      [x[i,j],y[i,j],z[i,j]])

    ax.plot_wireframe(x, y, z,  rstride=4, cstride=4, color='#2980b9', alpha=0.2)


def read_imaris_data(file_pattern):
    #Load numpy arrays
    co_position = np.load(file_pattern + '_crossover.npy')
    ids = np.load(file_pattern + '_ids.npy')
    positions = np.load(file_pattern + '_positions.npy')

    return co_position, ids, positions

def plot_ellipsoid(experiment_row):
    fig = plt.figure()

    ax = fig.add_subplot(111, projection='3d')
    ax.set_aspect("equal")

    # for fil in experiment_row['Positions']:
    #     fil = np.array(fil)
    #     ells = nestle.bounding_ellipsoids(fil)
    #     ax.plot(fil.T[0], fil.T[1], fil.T[2], "ko")
    #
    #     for ell in ells:
    #         plot_ellipsoid_3d(ell, ax)

    fil = experiment_row['Positions'][10]
    fil = np.array(fil)
    ells = nestle.bounding_ellipsoids(fil,0.0215)
    ax.plot(fil.T[0], fil.T[1], fil.T[2], "ko")

    for ell in ells:
        plot_ellipsoid_3d(ell, ax)

    plt.show()
    ax.set_xlabel("X")
    ax.set_xlim(ax.get_xlim()[::-1])
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    plt.show()

def plot_convex_hull(experiment_row):
    fig = plt.figure()

    ax = fig.add_subplot(111, projection='3d')
    ax.set_aspect("equal")

    for fil in experiment_row['Positions']:
        fil = np.array(fil)
        hull = ConvexHull(fil)
        ax.plot(fil.T[0], fil.T[1], fil.T[2], "ko")

        # 12 = 2 * 6 faces are the simplices (2 simplices per square face)
        for s in hull.simplices:
            s = np.append(s, s[0])  # Here we cycle back to the first coordinate
            ax.plot(fil[s, 0], fil[s, 1], fil[s, 2], "r-", alpha=0.2)


    plt.show()
    ax.set_xlabel("X")
    ax.set_xlim(ax.get_xlim()[::-1])
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    plt.show()

def add_volume_and_hull(df):
    positions = np.array(df['Positions'])
    convex_hull = []
    volume = []
    for elmnt in positions:
        poly = ConvexHull(elmnt)
        volume.append(poly.volume)
        convex_hull.append(poly)

    df['Hull'] = convex_hull
    df['Volume'] = volume

    return df
#########################
#Import and tidy experiment data
###########################
def load_data():
    exp_data = pd.read_csv('sample_info.csv')

    file_path = 'Jaso_data/'

    crossover = []
    ids = []
    positions = []

    for file_name in exp_data['Filename']:
        test_co, test_ids, test_positions = read_imaris_data(file_path + file_name)
        crossover.append(test_co)
        ids.append(test_ids)
        positions.append(test_positions)

    exp_data['Crossover'] = crossover
    exp_data['Ids'] = ids
    exp_data['Positions'] = positions

    lens = [len(item) for item in exp_data['Crossover']]
    df_out = pd.DataFrame( {"Filename" : np.repeat(exp_data['Filename'].values,lens),
                    "Stage" : np.repeat(exp_data['Stage'].values,lens),
                    "Genotype" : np.repeat(exp_data['Genotype'].values,lens),
                   "Crossover" : np.hstack(exp_data['Crossover']),
                   "Ids" : np.hstack(exp_data['Ids']),
                   "Positions" : np.hstack(exp_data['Positions'])
                  })

    chr_order = ['chr_X','chr_III','chr_V']

    df_out.sort_values(['Filename', 'Ids'])
    df_out['Chromosome'] = np.repeat(chr_order, len(df_out)/3)
    return df_out
