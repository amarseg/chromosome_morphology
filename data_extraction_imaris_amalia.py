###############

path='C:\\Users\\atamame\\Desktop\\Results_Imaris'

import ImarisLib
import numpy as np
import os
 

vLib = ImarisLib.ImarisLib()

vImaris = vLib.GetApplication(0)# Can run multiple Imaris applications – (0) denotes ID of Imaris application.

vImaris.GetVersion()

#Filament identification

vFilaments = vImaris.GetSurpassSelection();
vLabels = vFilaments.GetLabels;
filament_name = vFilaments.GetName()

vFilaments = vImaris.GetFactory().ToFilaments(vFilaments);

vFilamentNo=vFilaments.GetNumberOfFilaments();        

file_name = os.path.basename(vImaris.GetCurrentFileName())

positions = [];
dendrite_index =[];

for i in range(vFilamentNo):
	vPositionsXYZ = vFilaments.GetPositionsXYZ(i);  #XYZ coordinates of vertices for specified filament (0 – first filament in this case)
	positions.append(vPositionsXYZ);
	
	vBegin = vFilaments.GetBeginningVertexIndex(i);
	dendrite_index.append(vBegin);

pos_array = np.array(positions);
dendrite_array = np.array(dendrite_index);
ids_array = np.array(vFilaments.GetIds())
np.save(path +'\\' + file_name.split('.')[0] + filament_name + '_positions.npy',pos_array);
np.save(path +'\\' + file_name.split('.')[0] + filament_name + '_crossover.npy',dendrite_array);
np.save(path + '\\' +file_name.split('.')[0] + filament_name + '_ids.npy', ids_array);
################
