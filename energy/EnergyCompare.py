
def plot_energy_compare(ax,times,tcmin,label,names,energy_quant,energy_mult,colors,mean_to=10,\
						identify=[False],identify_ind=[None],identify_markers=[None]):
	"""
		In:
			ax : axes which to plot on (looped over each panel)
			times : time array ( len(times) == read_pkl(energy_quant[i]) )
			tcmin : min species cyclotron period 
			label : label to identify each panel (i.e. He3 concentration)
			energy_quant : energy quantity to load (KE or field)
			energy_mult : multipliers used for scaling (so consistent units J/m^3)
			colors : array of colors which to plot each quant
			identify : flag for whether identiftying times (for gyro-resonance)
			identify_ind : array of indices in time and dEnergy with which to set a scatter point (for gyro-resonance)
			indentify_markers : array of marker shapes which to scatter plot over data (for gyro-resonance)
		Out:
			ax : the same axes (panel) plotted on
	"""
	# if not identifying, make false array len of No. energies plotting
	if len(identify)==1:
		identify = [False]*len(energy_quant)
	# loop through each field and ion species
	for i in range(len(energy_quant)):
		Energy = read_pkl(energy_quant[i])*energy_mult[i]
		meanEnergy = np.mean(Energy[:mean_to])
		dEnergy = Energy-meanEnergy
		# print(len(dEnergy),len(times))
		ax.plot(times/tcmin,dEnergy,color=colors[i],zorder=1)#,label=names[i])
		if identify[i]:
			if not identify_markers[0]:
				identify_markers = ['x']*len(identify_ind)
			for j in range(len(identify_ind)):
				ax.scatter(times[identify_ind[j]]/tcmin,dEnergy[identify_ind[j]],marker=identify_markers[j],zorder=2,\
							s=50,label='_nolegend_',edgecolor='k',facecolor='none')#colors[i])
				ax.axvline(times[identify_ind[j]]/tcmin,color='k',linestyle='--',alpha=0.1,zorder=1,label='_nolegend_')
		ax.annotate(label,xy=(0.1,0.9),xycoords='axes fraction',ha='left',va='bottom')
	return ax

def energy_compare(sims,labels,tmax=7,colors=None,mean_to=10,frac=1,figname='',\
					identify_mat=[None],identify_indmat=[[None]],identify_markersmat=[[None]]):
	"""
	Plot a NxM panelled image of the energies of each sim in sims, with one label for each (if shared species)
		In:
			sims	: list of simulations to compare in an NxM plot
			labels	: labels annotated at the top left of each panel (i.e. He3 concentration)
			mean_to	: number of files (through time) to find mean initial energy of to find difference to
			tmax	: maximum time range to plot data (normalised to minority species from getIonSpecies(d0)[-1])
			frac  	: fraction of files to plot for the energies 
			figname : identifier to save fig(s)
			identify : array of booleans determining if to identify each energy quantity (i.e. if particle,particle,field then T,T,F)
			identify_indmat : 2d matrix of indices with which to identify the given data point (gyro-resonance) 
			identify_markersmat : 2d matrix of markers to use when identifying each point (gyro-resonance)
		Out:
			None
	"""
	# try 3, then 2 then 1x1 (i.e. single sim). If none of the above then opt for the least difference of 2 or 3
	try:
		M = N = 1
		n = np.array([3,2])
		for ni in n:
			m = int(len(sims))/ni
			print(m)
			if m%1==0:
				M = m ; N = ni
				break
			else: 
				M = 1 # single column
	except:
		N = M = 1 # single sim/image

	# check if self-defined or whether to use rainbow
	if not colors:
		colors = plt.cm.rainbow(np.linspace(0,1,5)) # two fields, three species

	# check if not identifying
	if identify_mat[0] == None:
		identify_mat	= [False for i in range(len(colors))]
		identify_indmat = [[None] for i in range(len(sims))]
		identify_markersmat = [[None] for i in range(len(sims))]
	if not identify_markersmat[0][0]:
		identify_markersmat = [[None] for i in range(len(sims))]

	# setup figure with N rows and M columns
	fig,axs=plt.subplots(nrows=int(M),ncols=int(N),figsize=(6,10),sharex=True,sharey=True)
	if N != 1 or M != 1: # multiple sims
		fig.subplots_adjust(hspace=0.075,wspace=0.075)
		axs=axs.ravel() # unravel axes
	else: print('# SINGLE SIMULATION #')

	# loop through each sim
	c=0
	for sim in sims:
		print(sim)
		# load sim & times
		simloc = getSimulation(sim)
		#_=energies(simloc)
		times = read_pkl('times')
		d0 = sdfread(0)

		# load field and particle names
		fields = ['Magnetic_Field_Bz','Electric_Field_Ex'] # hard-coded field plot
		efields = ['Bzenergy','Exenergy']
		ionspecies = list(filter(None,list(getIonSpecies(d0)))) # get all ion species (assuming all have same No. and type of species)
		tcmin = 2*const.PI/getCyclotronFreq(d0,ionspecies[-1]) # minority species
		energy_quant = efields + [i+'_KEdens' for i in ionspecies]

		# energy multipliers (for fields only) 
		energy_mult=[] ; names=[] ; ionlabels=[]
		for f in fields:
			_,l,m=Endict(f)
			energy_mult.append(m)
			names.append(l)
		for i in ionspecies:
			ionlabels.append(getIonlabel(i))		

		energy_mult = energy_mult + list(np.ones(len(ionspecies)))
		names = names + ionlabels # append names with labels of ions
		print(identify_indmat,identify_markersmat)
		print(identify_mat,identify_indmat[c],identify_markersmat[c])
		if N!=1 or M!=1: # multiple sims plot
			axs[c] = plot_energy_compare(axs[c],times,tcmin,labels[c],names,energy_quant,energy_mult,colors,identify=identify_mat,\
										identify_ind=identify_indmat[c],identify_markers=identify_markersmat[c])
			axs[c].axhline(0,linestyle='--',color='darkgrey')
		else: # single sim plot
			axs = plot_energy_compare(axs,times,tcmin,labels[0],names,energy_quant,energy_mult,colors)
			axs.axhline(0,linestyle='--',color='darkgrey')
		c+=1
		os.chdir('..')

	# formatting axes and limits
	if N!=1 or M!=1: # multiple sims
		axs[0].set_xlim(0,tmax)
		axs[0].set_ylim(-250,100)
		axs[0].locator_params(axis='y',nbins=4)
		# axs[len(axs)//2].set_ylabel(r'$\Delta u$'+'  ['+r'$Jm^{-3}$'+']',**tnrfont)
		# axs[0].set_ylabel(r'$\Delta u$'+'  ['+r'$Jm^{-3}$'+']',**tnrfont)
		# for n in range(int((M-1)*(N-1)+1),int(M*N)):
		# 	axs[n].set_xlabel(r'$t$'+getOmegaLabel(ionspecies[-1])+r'$/2\pi$',**tnrfont)
	else: # single sim
		axs.set_xlim(0,tmax)	
		axs.set_ylim(-250,100)
		# axs.set_ylabel(r'$\Delta u$'+'  ['+r'$Jm^{-3}$'+']',**tnrfont)
		# axs.set_xlabel(r'$t$'+getOmegaLabel(ionspecies[-1])+r'$/2\pi$',**tnrfont)
	legend = fig.legend(names,loc='upper center',ncol=len(energy_quant),bbox_to_anchor=(0.5,0.97),borderpad=0.1,\
						columnspacing=0.5,handlelength=1.5) # ncol=len(energy_quant)
	fig.supylabel(r'$\Delta u$'+'  ['+r'$Jm^{-3}$'+']',**tnrfont,x=-0.05)
	fig.supxlabel(r'$t$'+getOmegaLabel(ionspecies[-1])+r'$/2\pi$',**tnrfont,y=0.03)
	# plt.show()
	fig.savefig('energy_compare_{}.png'.format(figname),bbox_inches='tight')
	return None


if __name__=='__main__':
	from func_load import *
	from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
	import copy
	import GyroResonance as gr

	# # D-T
	# os.chdir('/storage/space2/phrmsf/traceT')
	# sims = np.sort([i for i in os.listdir() if 'traceT' in i])
	# sims = np.append(sims,sims[0])[1:] # largest to smallest tritium concentration (naming convention)
	# labels = [i[-2:] for i in sims]
	# # just zero case
	# sims = ['traceT_D_100_T_00']
	# labels = ['0']
	# energy_compare(sims,labels,colors=['b','g','r','orange','m'],tmax=7,figname='zero')

	# D-He3
	home = '/storage/space2/phrmsf/lowres_D_He3/'
	os.chdir(home)
	sims = np.sort([i for i in os.listdir(home) if 'p_90' in i])[1:] # excluding zero
	labels = [int(i[2:4]) for i in sims]
	# sims = np.array([home+i for i in sims])
	# just zero case
	# sims = ['0_00_p_90']
	# labels = ['0']

	# linear times
	identify_ind = [np.linspace(100,12000,8,dtype=int) for i in range(len(sims))]
	identify_markers = [['o','s','v','^','<','>','X','D'] for i in range(len(sims))]
	# energy comparison
	energy_compare(sims,labels,colors=['b','g','r','orange','m'],tmax=10,identify_mat=[False,False,True,True,False],\
					identify_indmat=identify_ind,identify_markersmat=identify_markers,figname='time_scatter')
	# gyro-resonance at times specified, single panel
	gr.gyro_time_compare(home,sims,identify_indmat=identify_ind,identify_markersmat=identify_markers,\
						multipanel=False,figname='singlepanel',plot_du=False,labels=labels)
	# multi-panel
	# gr.gyro_time_compare(home,sims,identify_indmat=identify_ind,identify_markersmat=identify_markers,\
	# 					multipanel=True,figname='multi_panel')

	# single sim energy compare
	# sims = ['0_00_p_90']
	# labels = ['0']
	# energy_compare(sims,labels,colors=['b','g','r','orange','m'],tmax=10,figname='zero') # ,figname='time_scatter')
