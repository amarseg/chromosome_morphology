import ImarisLib

vLib = ImarisLib.ImarisLib()
vImaris = vLib.GetApplication(0)                # Can run multiple Imaris applications – (0) denotes ID of Imaris application.
vImaris.GetVersion()

#Some testing script
vImage = vImaris.GetDataSet()
vSizeZ = vImage.GetSizeZ()
print vSizeZ, vImage.GetType()

#Filament identification
vFilaments = vImaris.GetSurpassSelection();
vFilaments.GetName();
vFilaments = vImaris.GetFactory().ToFilaments(vFilaments);

vIds = vFilaments.GetIds()

vFilamentNo=vFilaments.GetNumberOfFilaments();

filament_positions=[]
vPositionsXYZ = vFilaments.GetPositionsXYZ();  #XYZ coordinates of vertices for specified filament (0 – first filament in this case)
