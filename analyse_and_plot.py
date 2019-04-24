import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import os
import pandas as pd
from scipy.spatial import ConvexHull
import seaborn as sns
import nestle
import math
import vg
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

def plot_convex_hull(experiment_df, exp_id, convex = False):
    fig = plt.figure()

    ax = fig.add_subplot(111, projection='3d')
    ax.set_aspect("equal")

    subset = experiment_df.loc[experiment_df['Filename'] == 'EP_MD_4_TIRF- Filtered_Channel Alignment']

    for fil in subset['Positions']:
        fil = np.array(fil)
        hull = ConvexHull(fil)
        ax.plot(fil.T[0], fil.T[1], fil.T[2], "ko")

        if convex:
            # 12 = 2 * 6 faces are the simplices (2 simplices per square face)
            for s in hull.simplices:
                s = np.append(s, s[0])  # Here we cycle back to the first coordinate
                ax.plot(fil[s, 0], fil[s, 1], fil[s, 2], "r-", alpha=0.2)

    ax.set_xlabel("X")
    ax.set_xlim(ax.get_xlim()[::-1])
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    plt.show()

def add_volume_and_hull(df):
    positions = np.array(df['Positions'])
    convex_hull = []
    volume = []
    grad1 = []
    grad2 = []
    grad3 = []
    for elmnt in positions:
        poly = ConvexHull(elmnt)
        volume.append(poly.volume)
        convex_hull.append(poly)
        g = np.gradient(elmnt)
        grad1.append(g[0])
        grad2.append(g[1])
        g_g = np.gradient(g)
        grad3.append(g_g[0])

    df['Hull'] = convex_hull
    df['Volume'] = volume
    df['Gradient_1'] = grad1
    df['Gradient_2'] = grad2
    df['Gradient_3'] = grad3

    return df

def calculate_distance(p1, p2):
    differ = np.array(p1) - np.array(p2)
    dist = math.sqrt(sum(differ)**2)
    return dist

def calculate_total_distance(df):
    positions = np.array(df['Positions'])
    distances = []
    for fil in positions:
        ind_dist = []
        for n in range(0,len(fil)-1):
            d = calculate_distance(fil[n], fil[n+1])
            ind_dist.append(d)
        total_distance = sum(ind_dist)
        distances.append(total_distance)
    return distances

def angle_calculation(df):
    positions = np.array(df['Positions'])
    angles = []
    mean_angle = []
    for fil in positions:
        todo = []
        for n in range(0,len(fil)-2):
            v0 = np.array(fil[n]) - np.array(fil[n+1])
            v1 = np.array(fil[n+1]) - np.array(fil[n+2])
            # angle = np.math.atan2(np.linalg.det([v0,v1]),np.dot(v0,v1))
            # angle_deg = np.degree(angle_deg)
            angle = vg.angle(v0,v1)
            todo.append(angle)
        angles.append(todo)
        avg = np.median(todo)
        mean_angle.append(avg)
    return angles, mean_angle

def calculate_length_origin(df):
    positions = np.array(df['Positions'])
    start_end_dist = []
    for fil in positions:
        start = np.array(fil[0])
        end = np.array(fil[-1])
        start_end = calculate_distance(start, end)
        start_end_dist.append(start_end)
    return start_end_dist

def add_x_axis(df):
    x_axis = []
    for index, row in df.iterrows():
        indexes = list(range(1,len(row['Angles'])))
        poses =  [x - row['Crossover'] for x in indexes]
        x_axis.append(poses)
    return x_axis

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
    df_out['Chromosome'] = chr_order * 221
    df_out.to_csv('py_output.csv')
    return df_out


#########################
#Plotting and analysing data
###########################
df_out = load_data()
t = add_volume_and_hull(df_out)
print('Added volume and convex hull polygon to dataset')

total_distance = calculate_total_distance(t)
t['Total_distance'] = total_distance
print('Total Distance calculated and added to dataset')

d_start_end = calculate_length_origin(t)
t['Start_end_distance'] = d_start_end
t['Compaction_ratio'] = t['Total_distance']/t['Start_end_distance']
print('Distance between start and end of filament, plus compaction ratio added')

angles, mean_angle = angle_calculation(t)
t['Angles'] = angles
t['mean_angle'] = mean_angle

x_axis = add_x_axis(t)
t['x_axis'] = x_axis

print('Added angles between points')

t.to_csv('processed_data.csv')
print('Dataset saved')
