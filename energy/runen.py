from list_new import *
import my_constants as const

kwargs={'interpolation':'nearest','origin':'lower','aspect':'auto'}

def plotenergy(ax,times,tnorm,energy_quant,energy_mult,name,color='k',frac=1,mean_to=10):
	Energy=read_pkl(energy_quant)
	mean_Energy=np.mean(Energy[:mean_to])
	energy_plot = (Energy-mean_Energy)*energy_mult
	ax.plot(times[::frac]/tnorm,energy_plot[::frac],label=name,color=color)
	return ax
	
simlst = ['D_He3_0_10_COLD','D_He3_0_10_min_p_0_15','D_He3_0_10_min_p_0_9','D_He3_min_He4_0_8','D_He3_min_He4_0_9']
#simlst = ['D_He3_min_He4_0_8']
fig,axs = plt.subplots(nrows=len(simlst),ncols=3,sharex='col',sharey='col',gridspec_kw={'width_ratios': [1,1,1]})
colors=['b','cyan','g','r','m','orange','k']
maxt = 0
for i in range(len(simlst)):
	print(simlst[i])
	sim_loc = getSimulation('/storage/space2/phrmsf/'+simlst[i])
	times = read_pkl('times')
	d0 = sdfread(0)
	species = getIonSpecies(d0)
	dx = getdxyz(d0) ; dt = getdt() 
	wcD = const.qe*getMeanField3D(d0,'Magnetic_Field_B')/getMass('Deuterons')
	LDe = getDebyeLength(d0,'Electrons')	
#	## energy
#	energy_quant, names, energy_mult, fieldquant = getEnergyLabels(d0,species)
#	if i == 0:
#		maxt = times[-1]
#	else:
#		if maxt > times[-1]:
#			maxt = times[-1]
#	for j in range(len(energy_quant)):
#		axs[i,0] = plotenergy(axs[i,0],times,1,energy_quant[j],energy_mult[j],names[j],colors[j])
#	## FT2d
#	fm = load_batch_fieldmatrix([],'Magnetic_Field_Bz')
#	dfm = fm - np.mean(fm[0:10,:])
#	try:
#		FT2d = read_pkl('FT_2d_Magnetic_Field_Bz')
#	except:
#		FT2d = get2dTransform(dfm)
#		dumpfiles(FT2d,'FT_2d_Magnetic_Field_Bz')
#	print(wcD,LDe)
#	klim = 0.5*2*const.PI/dx
#	wlim = 0.5*2*const.PI/dt
#	axs[i,1].imshow(np.log10(FT2d),**kwargs,extent=[0,klim*LDe,0,wlim/wcD],cmap='magma',clim=(-2,6))
#	## bicoherence
	karea = 0.06
	warea = 40
	bicname = 'Bicohmat_ka_wa_{}_{}'.format(karea,warea)
	bisname = 'Bispecmat_ka_wa_{}_{}'.format(karea,warea)
	try:
		bic = read_pkl(bicname)
	except:
		fm = load_batch_fieldmatrix([],'Magnetic_Field_Bz')
		dfm = fm - np.mean(fm[0:10,:])
		## area
		verts = [
		   (0., 0.),  # left, bottom
		   (0.,karea),  # left, top
		   (warea, karea),  # right, top
		   (warea, 0.),  # right, bottom
		   (0., 0.),  # ignored
		]
		codes = [
			Path.MOVETO,
			Path.LINETO,
			Path.LINETO,
			Path.LINETO,
			Path.CLOSEPOLY,
		]
		area = Path(verts, codes)
		bis,bic=bispectrum2D(dfm,dt,getGridlen(d0),times[-1],area,len(times)//5,noverlap=len(times)//10,norm=[1/LDe,wcD],window=True,bispectrum=True)
		dumpfiles(bis,bisname)
		dumpfiles(bic,bicname)
	fig2,ax2,_,_=plot_bicoh(bic,extent=[0,karea,0,karea],bispectrum=False,smooth=True,cbar=True,clim=(0, 1),cmap ='jet')
	axs[i,2] = ax2


axs[0,0].set_xlim(0,maxt)
axs[0,1].set_xlim(0,0.075)
axs[0,1].set_ylim(0,50)

#for i in range(axs.shape[0]):
plt.show()
