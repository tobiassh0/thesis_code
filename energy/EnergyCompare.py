
def plot_energy_compare(ax,times,tcmin,label,names,energy_quant,energy_mult,colors,mean_to=10):
	# loop through each field and ion species
	for i in range(len(energy_quant)):
		Energy = read_pkl(energy_quant[i])*energy_mult[i]
		meanEnergy = np.mean(Energy[:mean_to])
		dEnergy = Energy-meanEnergy
		# print(len(dEnergy),len(times))
		ax.plot(times/tcmin,dEnergy,color=colors[i])#,label=names[i])
	ax.annotate(int(label),xy=(0.1,0.9),xycoords='axes fraction',ha='left',va='bottom')
	return ax

def energy_compare(sims,labels,tmax=7,colors=None,mean_to=10,frac=1,figname=''):
	"""
	Plot a NxM panelled image of the energies of each sim in sims, with one label for each (if shared species)
		In:
			sims		: 
			labels	: 
			mean_to	: 
			tmax		: 
			frac  	: fraction of files to plot for the energies 
		Out:
			
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

	# setup figure with N rows and M columns
	fig,axs=plt.subplots(nrows=int(N),ncols=int(M),figsize=(12,6),sharex=True,sharey=True)
	if N != 1 or M != 1: # multiple sims
		fig.subplots_adjust(hspace=0.075,wspace=0.075)
		axs=axs.ravel() # unravel axes
	else: print('# SINGLE SIMULATION #')

	# loop through each sim
	c=0
	for sim in sims:
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
		names = names + ionlabels
		if N!=1 or M!=1: # multiple sims plot
			axs[c] = plot_energy_compare(axs[c],times,tcmin,label[c],names,energy_quant,energy_mult,colors)
		else:# single sim plot
			axs = plot_energy_compare(axs,times,tcmin,labels[0],names,energy_quant,energy_mult,colors)
		# # loop through each field and ion species
		# for i in range(len(energy_quant)):
		# 	Energy = read_pkl(energy_quant[i])*energy_mult[i]
		# 	meanEnergy = np.mean(Energy[:mean_to])
		# 	dEnergy = Energy-meanEnergy
		# 	print(len(dEnergy),len(times))
		# 	axs[c].plot(times/tcmin,dEnergy,color=colors[i])#,label=names[i])
		# axs[c].annotate(int(labels[c]),xy=(0.1,0.9),xycoords='axes fraction',ha='left',va='bottom')
		c+=1

		os.chdir('..')

	# formatting axes and limits
	if N!=1 or M!=1: # multiple sims
		axs[0].set_xlim(0,tmax)
		axs[0].set_ylim(-250,100)
		axs[0].set_ylabel(r'$\Delta u$'+'  ['+r'$Jm^{-3}$'+']',**tnrfont)
		axs[len(axs)//2].set_ylabel(r'$\Delta u$'+'  ['+r'$Jm^{-3}$'+']',**tnrfont)
		for n in range(int((M-1)*(N-1)+1),int(M*N)):
			axs[n].set_xlabel(r'$t$'+getOmegaLabel(ionspecies[-1])+r'$/2\pi$',**tnrfont)
	else: # single sim
		axs.set_xlim(0,tmax)	
		axs.set_ylim(-250,100)
		axs.set_ylabel(r'$\Delta u$'+'  ['+r'$Jm^{-3}$'+']',**tnrfont)
		axs.set_xlabel(r'$t$'+getOmegaLabel(ionspecies[-1])+r'$/2\pi$',**tnrfont)
	legend = fig.legend(names,loc='upper center',ncol=len(energy_quant),bbox_to_anchor=(0.5,1.0),borderpad=0.1)
	fig.savefig('energy_compare_{}.png'.format(figname),bbox_inches='tight')
	return None

if __name__=='__main__':
	from func_load import *
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
	os.chdir('/storage/space2/phrmsf/lowres_D_He3')
	sims = np.sort([i for i in os.listdir() if 'p_90' in i])[1:] # excluding zero
	labels = [i[2:4] for i in sims]
	# just zero case
	sims = ['0_00_p_90']
	labels = ['0']
	energy_compare(sims,labels,colors=['b','g','r','orange','m'],tmax=10,figname='zero')
