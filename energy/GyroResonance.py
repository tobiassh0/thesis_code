
from func_load import *
import numpy as np
import matplotlib.pyplot as plt

def PlotGyroResonance(home,duarr,sims,species=['Deuterons','Tritons'],norm_spec='Alphas',through_time=False,identify=False,lims=((0,1),(0,1)),labels='',\
					  ylabel=r'$[\Delta u_1/\Delta u_2]_{max}$',xlabel=r'$(\xi_1/\xi_2)(m_2/m_1)(q_1/q_2)^2$',time_norm=r'$\tau_{cD}$'):
	"""
		IN:
			duarr 			: shape (len(sims), len(species), len(times))
			ylabel			: self exp.
			xlabel			: self exp.
			labels			: the labels of each simulation (i.e. concentration %), valid if identify=True
			identify		: identifies each simulation and annotates/plots a legend box (default False)
			through_time	: determines whether you plot a through_time or 1:1 du ratio (default False)
		OUT:
			None
	"""
	if identify:
		colors = cm.rainbow(np.linspace(0,1,len(sims)))
	else:
		colors = ['b']*len(sims)
		labels =  ''
	if not through_time:
		fig,ax=plt.subplots(figsize=(6,4))
		ax.plot(lims[0],lims[0],color='darkgray',linestyle='--') # 1:1 line
	else:
		fig,ax=plt.subplots(figsize=(8,6))

	_=getSimulation(home+sims[0])
	d0 = sdfread(0)
	tcnorm = 2*const.PI/getCyclotronFreq(d0,norm_spec)
	# mass and charge should be constant throughout sims (same species being analysed)
	marr = [getMass(sp) for sp in species]
	qarr = [getChargeNum(sp) for sp in species]
	c=0
	for j in range(len(sims)):
		_=getSimulation(home+sims[j])
		d0 = sdfread(0)
		times = read_pkl('times')
		xi1,xi2,_=getConcentrationRatios(d0)
		print('Concentration xi2 :',xi2)
		# plotting through time or statically at maximal change
		if not through_time:
			# maximal du per species 0 & 1
			maxdu0 = np.max(duarr[j,0,:])
			maxdu1 = np.max(duarr[j,1,:])
			ax.scatter((xi1/xi2)*(marr[1]/marr[0])*((qarr[0]/qarr[1])**2),maxdu0/maxdu1,color=colors[c])
			xy = ((xi1/xi2)*(marr[1]/marr[0])*((qarr[0]/qarr[1])**2),maxdu0/maxdu1)
		else: # through time for each sim
			ax.plot(times/tcnorm,np.abs((xi2/xi1)*duarr[j,0,:]/duarr[j,1,:]),color=colors[c])
			ax.set_ylim(1e-3,1e3)
			ax.set_xlim(0,times[-1]/tcnorm)
			ax.set_yscale('log')
		c+=1
	# saving and formatting
	if not through_time:
		if identify:
			import matplotlib as mpl
			mpl.rc('font',family='Times New Roman')
			ax.legend(labels,loc='best',fontsize=16)
			name = home+'/du_peak_vs_theory_label.png'
		else:
			name = home+'/du_peak_vs_theory.png'
		ax.set_ylabel(ylabel,**tnrfont) ; ax.set_xlabel(xlabel,**tnrfont)
		ax.set_xlim(lims[0]) ; ax.set_ylim(lims[1])		
		# fig.savefig(name,bbox_inches='tight')
	else:
		if identify:
			ax.legend(labels,loc='best') #,ncol=len(labels))
		ax.axhline((marr[1]/marr[0])*(qarr[0]/qarr[1])**2,linestyle='--',color='k')
		ax.set_ylabel(ylabel,fontsize=24)
		ax.set_xlabel(r'$t/$'+time_norm,**tnrfont)
		fig.savefig(home+'/du_ratio_vs_time.png',bbox_inches='tight')
	plt.show()
	return None

# Plots and shows the experiment vs theory plot for the ratio between two species change in energy density
# Also has the ability to plot du per-particle ratio (per species) through time for each simulation (nrows) -- will need to un-comment these lines
def majIons_edens_ratio(home,sims,species=['Deuterons','Tritons'],norm_spec='Alphas',plot=False,mean_to=10,N=50):
	"""
		Function to find the ratio between energy densities of two species (specified) across a range of simulations (specified) and plot it 
		either via a 1:1 correlation (through_time=False) or through time for each simulation (through_time=True). Is able to identify each
		simulation for either scenario (identify=True/False). Then saves each figure in the home directory that the code is executed in.
	
		IN:
			sims	 : the simulations you want to loop through and extract their du [J m^-3]
			species	 : the species present in each simulation which you want to compare
			norm_spec: species to normalise time to
			lims	 : the limits of the 1:1 line plot (default is between 0 and 1)
			plot	 : boolean determining whether to plot or not (default False)
			N		 : Convolution smoothing
			mean_to  : No. time files to calculate mean energy (DeltaE = E - meanE)
		OUT:
			duarr	 : array of change in energy densities through time of 
	"""	
	# mass and charge should be constant throughout sims (same species being analysed)
	if species == []:
		species = getIonSpecies(d0)[:-1] # remove last minority particle
	marr = [getMass(sp) for sp in species]
	qarr = [getChargeNum(sp) for sp in species]
	sim0 = sims[0]
	_=getSimulation(home+sim0)
	times0 = read_pkl('times')
	duarr = np.zeros((len(sims),len(species),len(times0)))
	for j in range(len(sims)):
		sim_loc = getSimulation(home+sims[j])
		d0 = sdfread(0)
		tcnorm = 2*const.PI/getCyclotronFreq(d0,norm_spec)
		## field & species energies
		xi1,xi2,_=getConcentrationRatios(d0)
		xi2_xi1 = xi2/xi1
		print(sims[j])
		print('Secondary conc.:: ',xi2)
		for i in range(len(species)):
			Energypart = read_pkl(species[i]+'_KEdens')/(1000*const.qe) # keV/m^3
			meanEnergypart = np.mean(Energypart[:mean_to])
			# Energypart = np.convolve(Energypart,np.ones(N)/N,mode='valid')
			du = (Energypart-meanEnergypart)
			duarr[j,i,:] = du
	if plot:
		PlotGyroResonance(home,duarr,sims,species,norm_spec=norm_spec)
	
	return duarr

if __name__=='__main__':
	from func_load import *	
	import matplotlib.cm as cm
	
	# # D-T
	# home = '/storage/space2/phrmsf/traceT/'
	# sims = np.sort([i for i in os.listdir(home) if 'traceT' in i])
	# sims = sims[1:]#np.append(sims,sims[0])[1:] # move 0% tritium to end
	# hlabels = np.array([int(i[-2:]) for i in sims])
	# # standard
	# ylabel=r'$[\Delta u_{T}/\Delta u_{D}]_{max}$'
	# xlabel=r'$(\xi_{T}/\xi_{D})(m_{D}/m_{T})(q_{T}/q_{D})^2$'
	# # # through time
	# # ylabel=r'$\left(\frac{\xi_D}{\xi_T}\right)\left|\frac{\Delta u_T(t)}{\Delta u_D(t)}\right|$'
	# duarr = majIons_edens_ratio(home,sims,species=['Deuterons','Tritons'],norm_spec='Alphas')
	# print(duarr,duarr.shape)
	# PlotGyroResonance(home,duarr,sims,species=['Deuterons','Tritons'],norm_spec='Alphas',labels=hlabels,\
	# 				  ylabel=ylabel,xlabel=xlabel,lims=((0,1),(0,1)),through_time=False,identify=False,time_norm=r'$\tau_{cD}$')
	
	# D-He3
	home = '/storage/space2/phrmsf/lowres_D_He3/'
	sims = np.sort([i for i in os.listdir(home) if 'p_90' in i])[1:]
	print(sims)
	hlabels = np.array([int(i[2:4]) for i in sims])
	# standard
	xlabel=r'$(\xi_{D}/\xi_{He3})(m_{He3}/m_{D})(q_{D}/q_{He3})^2$'
	ylabel=r'$[\Delta u_{D}/\Delta u_{He3}]_{max}$'
	# # through-time
	# ylabel=r'$(\Delta u_{D}/\Delta u_{He3})$' + ' ' + r'$(\xi_{He3}/\xi_D)$'
	duarr = majIons_edens_ratio(home,sims,species=['Deuterons','He3'],norm_spec='Protons')
	print(duarr,duarr.shape)
	PlotGyroResonance(home,duarr,sims,species=['Deuterons','He3'],norm_spec='Protons',labels=hlabels,\
					  ylabel=ylabel,xlabel=xlabel,lims=((0,8),(0,8)),through_time=True,identify=True,time_norm=r'$\tau_{cp}$')
