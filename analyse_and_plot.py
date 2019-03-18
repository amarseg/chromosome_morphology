import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D


def read_imaris_data(file_pattern):
    #Load numpy arrays
    co_position = np.load(file_pattern + '_crossover.npy')
    ids = np.load(file_pattern + '_ids.npy')
    positions = np.load(file_pattern + '_positions.npy')

    return co_position, ids, positions

file_test = 'Jaso_data/EP_MD_4_TIRF- Filtered_Channel Alignment'
test_co, test_ids, test_positions = read_imaris_data(file_test)


print(test_ids)
print(test_co)

fig = plt.figure()

ax = fig.add_subplot(111, projection='3d')

for fil in test_positions:
    t = np.array(fil)
    print(t.gradient())
    ax.plot(t[:,0],t[:,1],t[:,2])

plt.show()
