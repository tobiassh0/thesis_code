'''
	Example script to highlight usage of functions for EPOCH analysis.
	
	Here we load the packages through func_load, get the simulation location 
	as a path and change to that dir. Then we list all sdf files as a t_files
	array. The quantity E_x electric field is then loaded linearly through time 
	and space into the array Exfield. This is then plotted as an un-normalised
	heat-map (norm=1 default) and returned to this script so it can be saved. 
'''

from func_load import *
# loads all functions you might need

# get simulation location
simLocation = getSimulation('/sim/file/location') # so that the dir can read each .sdf file

# list of all output files
t_files = list_sdf(simLocation)

# load fieldmatrix of the field 'Ex' (loads in linear not parallel)
quantity = 'Electric_Field_Ex'
Exfield = load_batch_fieldmatrix(t_files,quantity,para=False)

# get simulation length	and duration 
file0 =	sdfread(0)
L = getGridlen(file0)
times = read_pkl('times')
T = times[-1]

# plot field
xlabel='Position, '+'$m$'
ylabel='Time, '+r'$t$'
cbar_label=r'$E_x$'+'  '+r'$[V/m]$'
xlim = L ; ylim = T
fig,ax=plotNormHeatmap(matrix,xlabel=xlabel,ylabel=ylabel,cbar_label=\
		cbar_label,xylim=((0,L),(0,T)),norm=1)

# save image
fig.savefig('EX-heatmap.png')


