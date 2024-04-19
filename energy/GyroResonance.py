
from func_load import *

# Plots and shows the experiment vs theory plot for the ratio between two species change in energy density
# Also has the ability to plot du per-particle ratio (per species) through time for each simulation (nrows) -- will need to un-comment these lines
def majIons_edens_ratio(sims,species=['Deuterons','Tritons'],time_norm=r'$\tau_{cD}$',ylabel=r'$[\Delta u_1/\Delta u_2]_{max}$',\
								xlabel=r'$(\xi_1/\xi_2)(m_2/m_1)(q_1/q_2)^2$',tmax=6.1,labels=[None],lims=((0,1),(0,1)),\
								identify=False,through_time=False,plot=True):
	"""
		Function to find the ratio between energy densities of two species (specified) across a range of simulations (specified) and plot it 
		either via a 1:1 correlation (through_time=False) or through time for each simulation (through_time=True). Is able to identify each
		simulation for either scenario (identify=True/False). Then saves each figure in the home directory that the code is executed in.
		
		params in:
			sims				: the simulations you want to loop through and extract their du [J m^-3]
			species			: the species present in each simulation which you want to compare
			time_norm		: normalised time with which to plot the through_time version
			ylabel			: self exp.
			xlabel			: self exp.
			labels			: the labels of each simulation (i.e. concentration %), valid if identify=True
			tmax 				: the maximal time to find the maximal du, in units of time_norm 
			lims				: the limits of the 1:1 line plot (default is between 0 and 1)
			identify			: identifies each simulation and annotates/plots a legend box (default False)
			through_time	: determines whether you plot a through_time or 1:1 du ratio (default False)
			plot				: boolean determining whether the 
		params out:
			None, plots the respective figures 
	"""	
	mean_to = 10
	N=50
	c=0
	if plot:
		if identify:
			colors = cm.rainbow(np.linspace(0,1,len(sims)))
		else:
			colors = ['b']*len(sims)
			# labels = [None]*len(sims)
		if not through_time:
			fig,ax=plt.subplots(figsize=(6,4))
			ax.plot(lims[0],lims[0],color='darkgray',linestyle='--') # 1:1 line
		else:
			fig,ax=plt.subplots(figsize=(6,6))	
	home = os.getcwd()
	xiarr=np.zeros((len(sims),len(species)))
	# mass and charge should be constant throughout sims (same species being analysed)
	marr = [getMass(sp) for sp in species]
	qarr = [getChargeNum(sp) for sp in species]
	for j in range(len(sims)):
		sim_loc = getSimulation(sims[j])
		times = read_pkl('times')
		d0 = sdfread(0)
		nx = len(getQuantity1d(d0,'Derived_Number_Density'))
		n0 = getQuantity1d(d0,'Derived_Number_Density_Electrons')
		wc1 = getCyclotronFreq(d0,species[0])
		tc1 = 2*const.PI/wc1
		dt = (times[-1]-times[0])/len(times)
		dt_p = dt/tc1
		## field & species energies
		maxdu = np.zeros(len(species)) # reset per sim
		xi1,xi2,_=getConcentrationRatios(d0)
		xi2_xi1 = xi2/xi1
		xiarr[j,:]=np.array([xi1,xi2])
		duarr=[]
		print(sims[j])
		print('Secondary conc.:: ',xi2)
		for i in range(len(species)):
			Energypart = read_pkl(species[i]+'_KEdens')/(1000*const.qe) # keV/m^3
			meanEnergypart = np.mean(Energypart[:mean_to])
			Energypart = np.convolve(Energypart,np.ones(N)/N,mode='valid')
			timespart = np.linspace(0,max(times),len(Energypart))
			thresh = timespart/tc1 < tmax
			timespart = timespart[thresh]
			du = (Energypart[thresh]-meanEnergypart)
			maxdu[i] = np.max(du)
			duarr.append(du)
		duarr = np.array(duarr) # size (species, thresh time)
		if plot:
			if not through_time:
				ax.scatter((xi1/xi2)*(marr[1]/marr[0])*((qarr[0]/qarr[1])**2),maxdu[0]/maxdu[1],color=colors[c],label=labels[c])
				xy = ((xi1/xi2)*(marr[1]/marr[0])*((qarr[0]/qarr[1])**2),maxdu[0]/maxdu[1])
				# ax.annotate(str(labels[c]),xy=xy,xycoords='data',xytext=(.5,0),textcoords='offset fontsize',fontsize=16,fontname='Times New Roman')
			else:	# through time for each sim
				ax.plot(timespart/tc1,np.abs((1/xi2_xi1)*duarr[1]/duarr[0]),color=colors[c],label=str(labels[c])+'%')
				# ax[c].annotate(str(labels[c])+'%',xy=(0.9,0.85),xycoords='axes fraction',**tnrfont)
				ax.set_ylim(1e-3,1e3)
				ax.set_xlim(0,tmax)
				ax.axhline(marr[0]/marr[1],linestyle='--',color='k')
				ax.set_yscale('log')
		os.chdir(home)
		c+=1
	if plot:
		if through_time:
			if identify:
				plt.legend(loc='best')#,ncol=len(labels))
			ax.set_ylabel(ylabel,fontsize=24)
			ax.set_xlabel(xlabel+time_norm,**tnrfont)
			fig.savefig('du_ratio_vs_time.png',bbox_inches='tight')
		else:
			if identify:
				import matplotlib as mpl
				mpl.rc('font',family='Times New Roman')
				ax.legend(loc='best',fontsize=16)
				name = 'du_peak_vs_theory_label.png'
			else:
				name = 'du_peak_vs_theory.png'
			ax.set_ylabel(ylabel,**tnrfont) ; ax.set_xlabel(xlabel,**tnrfont)
			ax.set_xlim(lims[0]) ; ax.set_ylim(lims[1])
			fig.savefig(name,bbox_inches='tight')
			plt.show()
	return duarr, xiarr, marr, qarr

if __name__=='__main__':
	from func_load import *	
	import matplotlib.cm as cm
	
	# # D-T
	# os.chdir('/storage/space2/phrmsf/traceT/')
	# sims = np.sort([i for i in os.listdir() if 'traceT' in i])
	# sims = sims[1:]#np.append(sims,sims[0])[1:] # move 0% tritium to end
	# hlabels = np.array([int(i[-2:]) for i in sims])
	# # # max du 
	# # ylabel=r'$[\Delta u_{T}/\Delta u_{D}]_{max}$'
	# # xlabel=r'$(\xi_{T}/\xi_{D})(m_{D}/m_{T})(q_{T}/q_{D})^2$'
	# # through time
	# ylabel=r'$\left(\frac{\xi_D}{\xi_T}\right)\left|\frac{\Delta u_T(t)}{\Delta u_D(t)}\right|$'
	# xlabel=r'Time,  '
	# _,_,_,_=majIons_edens_ratio(sims,species=['Deuterons','Tritons'],time_norm=r'$\tau_{cD}$',labels=hlabels,\
	# 									 xlabel=xlabel,ylabel=ylabel,lims=((0,1),(0,1)),identify=True,plot=True,through_time=True)
	
	# D-He3
	os.chdir('/storage/space2/phrmsf/lowres_D_He3/')
	sims = np.sort([i for i in os.listdir() if 'p_90' in i])[1:]
	print(sims)
	hlabels = np.array([int(i[2:4]) for i in sims])
	ylabel=r'$[\Delta u_{D}/\Delta u_{He3}]_{max}$'
	xlabel=r'$(\xi_{D}/\xi_{He3})(m_{He3}/m_{D})(q_{D}/q_{He3})^2$'
	_,_,_,_=majIons_edens_ratio(sims,species=['Deuterons','He3'],time_norm=r'$\tau_{cD}$',labels=hlabels,\
										 xlabel=xlabel,ylabel=ylabel,lims=((0,10),(0,10)),identify=False,plot=True,through_time=True)