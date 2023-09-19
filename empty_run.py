from func_load import *
from matplotlib.gridspec import GridSpec


sim_lst = ['lowres_D_He3/0_25_p_90']
quant = 'Magnetic_Field_Bz'

for sim in sim_lst:
	loc = getSimulation('/storage/space2/phrmsf/lowres_D_He3/0_34_p_90')
	d0 = sdfread(0)
	times = read_pkl('times')
	# cig plot
	ciggies(loc,species_lst=['Protons'])
	plt.show()
