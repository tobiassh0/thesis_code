from func_load import *
from matplotlib.gridspec import GridSpec


sim_lst = ['traceT_0_00']

for sim in sim_lst:
	loc = getSimulation('/storage/space2/phrmsf/traceT_0_00')
	L = getGridlen(sdfread(0))
	FT1d = read_pkl('FT_1d_Magnetic_Field_Bz')
#	field = load_batch_fieldmatrix('Magnetic_Field_Bz')
#	field = field - np.mean(field[0:10,:])
	FT_1d = read_pkl('FT_1d_test') # get1dTransform(field)
#	dumpfiles(FT_1d,'FT_1d_test')
	print(FT1d, FT_1d)
	vA = getAlfvenVel(sdfread(0))
	times = read_pkl('times')
	T = times[-1]
	dx = getdxyz(sdfread(0))
	dt = T/len(times)
	wcD = getCyclotronFreq(sdfread(0),'Deuterons')
	tcD = 2*const.PI/wcD

klim = 0.5*2*const.PI/dx ; tlim = T
extent = [0,klim/(wcD/vA),0,tlim/tcD]
print(extent)
fig = plt.figure()
fig.add_subplot(211)
plt.imshow(np.log10(FT1d),extent=extent,clim=(-2.5,4.0),**kwargs,cmap='seismic')
plt.title('loaded')
plt.xlim(0,40) ; plt.ylim(0,7)
fig.add_subplot(212)
plt.imshow(np.log10(FT_1d),extent=extent,clim=(-2.5,4.0),**kwargs,cmap='seismic')
plt.xlim(0,40) ; plt.ylim(0,7)
plt.title('calc dBz')
plt.savefig('FT_1d_comparison.png')


