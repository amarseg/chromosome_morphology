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

###########################
#Animation code (https://zulko.wordpress.com/2012/09/29/animate-your-3d-plots-with-pythons-matplotlib/)
############################
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import axes3d
import os, sys
import numpy as np


##### TO CREATE A SERIES OF PICTURES

def make_views(ax,angles,elevation=None, width=16, height = 12,
                prefix='tmprot_',**kwargs):
    """
    Makes jpeg pictures of the given 3d ax, with different angles.
    Args:
        ax (3D axis): te ax
        angles (list): the list of angles (in degree) under which to
                       take the picture.
        width,height (float): size, in inches, of the output images.
        prefix (str): prefix for the files created.

    Returns: the list of files created (for later removal)
    """

    files = []
    ax.figure.set_size_inches(width,height)

    for i,angle in enumerate(angles):

        ax.view_init(elev = elevation, azim=angle)
        fname = '%s%03d.jpeg'%(prefix,i)
        ax.figure.savefig(fname)
        files.append(fname)

    return files



##### TO TRANSFORM THE SERIES OF PICTURE INTO AN ANIMATION

def make_movie(files,output, fps=10,bitrate=1800,**kwargs):
    """
    Uses mencoder, produces a .mp4/.ogv/... movie from a list of
    picture files.
    """

    output_name, output_ext = os.path.splitext(output)
    command = { '.mp4' : 'mencoder "mf://%s" -mf fps=%d -o %s.mp4 -ovc lavc\
                         -lavcopts vcodec=msmpeg4v2:vbitrate=%d'
                         %(",".join(files),fps,output_name,bitrate)}

    command['.ogv'] = command['.mp4'] + '; ffmpeg -i %s.mp4 -r %d %s'%(output_name,fps,output)

    print(command[output_ext])
    output_ext = os.path.splitext(output)[1]
    os.system(command[output_ext])



def make_gif(files,output,delay=100, repeat=True,**kwargs):
    """
    Uses imageMagick to produce an animated .gif from a list of
    picture files.
    """

    loop = -1 if repeat else 0
    os.system('convert -delay %d -loop %d %s %s'
              %(delay,loop," ".join(files),output))




def make_strip(files,output,**kwargs):
    """
    Uses imageMagick to produce a .jpeg strip from a list of
    picture files.
    """

    os.system('montage -tile 1x -geometry +0+0 %s %s'%(" ".join(files),output))



##### MAIN FUNCTION

def rotanimate(ax, angles, output, **kwargs):
    """
    Produces an animation (.mp4,.ogv,.gif,.jpeg,.png) from a 3D plot on
    a 3D ax

    Args:
        ax (3D axis): the ax containing the plot of interest
        angles (list): the list of angles (in degree) under which to
                       show the plot.
        output : name of the output file. The extension determines the
                 kind of animation used.
        **kwargs:
            - width : in inches
            - heigth: in inches
            - framerate : frames per second
            - delay : delay between frames in milliseconds
            - repeat : True or False (.gif only)
    """

    output_ext = os.path.splitext(output)[1]

    files = make_views(ax,angles, **kwargs)

    D = { '.mp4' : make_movie,
          '.ogv' : make_movie,
          '.gif': make_gif ,
          '.jpeg': make_strip,
          '.png':make_strip}

    D[output_ext](files,output,**kwargs)

    for f in files:
        os.remove(f)
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
    scale_x = 2
    scale_y = 2
    scale_z = 2
    #ax.get_proj = lambda: np.dot(Axes3D.get_proj(ax), np.diag([scale_x, scale_y, scale_z, 1]))
    #ax.autoscale()
    subset = experiment_df.loc[experiment_df['Filename'] == 'EP_MD_4_TIRF- Filtered_Channel Alignment']

    fil = np.array(subset['Positions'][0])
    hull = ConvexHull(fil)
    ax.plot(fil.T[0], fil.T[1], fil.T[2], "ko")

    # for fil in subset['Positions']:
    #     fil = np.array(fil)
    #     hull = ConvexHull(fil)
    #     ax.plot(fil.T[0], fil.T[1], fil.T[2], "ko")
    #
    #     if convex:
    #         # 12 = 2 * 6 faces are the simplices (2 simplices per square face)
    #         for s in hull.simplices:
    #             s = np.append(s, s[0])  # Here we cycle back to the first coordinate
    #             ax.plot(fil[s, 0], fil[s, 1], fil[s, 2], "r-", alpha=0.2)

    ax.set_xlabel("X")
    ax.set_xlim(ax.get_xlim()[::-1])
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    #plt.show()
    angles = np.linspace(0,360,21)[:-1]
    rotanimate(ax, angles,'movie.gif',delay=20)

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

def filament_total_distance(fil):
    ind_dist = []
    for n in range(0,len(fil)-1):
        d = calculate_distance(fil[n], fil[n+1])
        ind_dist.append(d)
    total_distance = sum(ind_dist)
    return total_distance

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

def split_angles_into_short_and_long(df):
    positions = np.array(df['Angles'])
    co_pos = np.array(df['Crossover']) - 1
    short = []
    long = []
    for i in range(len(positions)):
        pos = positions[i]
        co = co_pos[i]
        arm1 = list(pos[0:co])
        arm2 = list(pos[co+1:len(pos)])
        l_arm1 = len(arm1)
        l_arm2 = len(arm2)

        if(l_arm1 > l_arm2):
            short.append(arm2)
            long.append(arm1)
        else:
            short.append(arm1)
            long.append(arm2)

    return short, long

def split_lengths_into_short_and_long(df):
    positions = np.array(df['Positions'])
    co_pos = np.array(df['Crossover'])
    short = []
    long = []
    for i in range(len(positions)):
        pos = positions[i]
        co = co_pos[i]
        arm1 = list(pos[0:co])
        arm2 = list(pos[co+1:len(pos)])
        l_arm1 = filament_total_distance(arm1)
        l_arm2 = filament_total_distance(arm2)

        if(l_arm1 > l_arm2):
            short.append(l_arm2)
            long.append(l_arm1)
        else:
            short.append(l_arm1)
            long.append(l_arm2)

    return short, long

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
t['sum_angle'] = list(map(sum, t['Angles']))/t['Total_distance']

x_axis = add_x_axis(t)
t['x_axis'] = x_axis

print('Added angles between points')

t.to_csv('processed_data.csv')
print('Dataset saved')

only_LP = t[t.Stage == 'LP']

short, long = split_angles_into_short_and_long(only_LP)
only_LP['Short_angle'] = short
only_LP['Long_angle'] = long
short_l, long_l = split_lengths_into_short_and_long(only_LP)
only_LP['Long_length'] = long_l
only_LP['Short_length'] = short_l

only_LP['Sum_angle_short'] = list(map(sum, only_LP['Short_angle']))/only_LP['Short_length']
only_LP['Sum_angle_long'] = list(map(sum, only_LP['Long_angle']))/only_LP['Long_length']
only_LP['Long_short_ratio'] = only_LP['Sum_angle_long']/only_LP['Sum_angle_short']
only_LP.to_csv('processed_LP.csv')


plot_convex_hull(t, 'EP_MD_4_TIRF- Filtered_Channel Alignment')
