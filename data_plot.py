import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import os
import pandas as pd
from scipy.spatial import ConvexHull
import seaborn as sns
import nestle
import math
#######
#Functions
######
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

#################
#Actual plotting
####################
t = pd.read_csv('processed_data.csv')
plt.figure(figsize=(16, 12))
sns.swarmplot(x="Genotype", y="Volume", palette=["r", "c", "y"], data=t)
plt.savefig('fig_3d/Swarm_Volume_genotype.png')

sns.boxplot(x='Chromosome', y='Volume', hue='Genotype',
            palette='Set3', data=t)
plt.savefig('fig_3d/Volume_genotype.png')

sns.catplot(x="Chromosome", y="Volume",
            hue="Genotype", col="Stage",
            data=t, kind="box",
            height=10, aspect=.7)
plt.savefig('fig_3d/Volume_stage.png')

sns.catplot(x="Chromosome", y="Volume",
            hue="Stage", col="Genotype",
            data=t, kind="box", palette = 'Set2',
            height=10, aspect=.7)
plt.savefig('fig_3d/facet_Volume_stage.png')

sns.catplot(x="Chromosome", y="Compaction_ratio",
            hue="Stage", col="Genotype",
            data=t, kind="box", palette = 'Set2',
            height=10, aspect=.7, showfliers=False)

plt.savefig('fig_3d/facet_compact_stage.png')

sns.catplot(x="Chromosome", y="Compaction_ratio",
            hue="Genotype", col="Stage",
            data=t, kind="box",
            height=10, aspect=.7, showfliers=False)

plt.savefig('fig_3d/facet_compact_genotype.png')

plt.show()
plot_convex_hull(t, 'EP_MD_4_TIRF- Filtered_Channel Alignment')
