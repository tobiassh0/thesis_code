
from func_load import *
import numpy as np
import matplotlib.pyplot as plt

def PlotGyroResonance(home,duarr,sims,species=['Deuterons','Tritons'],norm_spec='Alphas',through_time=False,identify=False,lims=((0,1),(0,1)),labels='',\
					  ylabel=r'$[\Delta u_1/\Delta u_2]_{max}$',xlabel=r'$(\xi_1/\xi_2)(m_2/m_1)(q_1/q_2)^2$',time_norm=r'$\tau_{cD}$',figname=''):
	"""
		IN:
			duarr 			: shape (len(sims), len(species), len(times))
			ylabel			: self exp.
			xlabel			: self exp.
			labels			: the labels of each simulation (i.e. concentration %), valid if identify=True
			identify		: boolean : identifies each simulation and annotates/plots a legend box (default False)
			through_time	: boolean : determines whether you plot a through_time or 1:1 du ratio (default False)
			lims 			: limits on the gyro-resonance theory plot, should be equal in x and y
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
			name = home+'/du_peak_vs_theory_label_{}.png'.format(figname)
		else:
			name = home+'/du_peak_vs_theory_{}.png'.format(figname)
		ax.set_ylabel(ylabel,**tnrfont) ; ax.set_xlabel(xlabel,**tnrfont)
		ax.set_xlim(lims[0]) ; ax.set_ylim(lims[1])		
		fig.savefig(name,bbox_inches='tight')
	else:
		if identify:
			ax.legend(labels,loc='best') #,ncol=len(labels))
		ax.axhline((marr[1]/marr[0])*(qarr[0]/qarr[1])**2,linestyle='--',color='k')
		ax.set_ylabel(ylabel,fontsize=24)
		ax.set_xlabel(r'$t/$'+time_norm,**tnrfont)
		fig.savefig(home+'/du_ratio_vs_time_{}.png'.format(figname),bbox_inches='tight')
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
			Energypart = read_pkl(species[i]+'_KEdens') # /(1000*const.qe) # J/m^3
			meanEnergypart = np.mean(Energypart[:mean_to])
			# Energypart = np.convolve(Energypart,np.ones(N)/N,mode='valid')
			du = (Energypart-meanEnergypart)
			duarr[j,i,:] = du
	if plot:
		PlotGyroResonance(home,duarr,sims,species,norm_spec=norm_spec)
	# return array of du for each sim, species and time
	return duarr

# plots the energy traces multiplied by the charge to mass gyro-ratio
def du_gyroratio(home,duarr,sims,species,norm_spec='Alphas',labels=''):
	nsims, nspec, nt = duarr.shape
	_=getSimulation(home+sims[0])
	tcnorm = 2*const.PI/getCyclotronFreq(sdfread(0),norm_spec)
	marr = [getMass(sp) for sp in species]
	qarr = [getChargeNum(sp) for sp in species]
	N = 2 ; M = 4
	# setup figure
	fig,ax=plt.subplots(figsize=(10,5),nrows=N,ncols=M,sharex=True,sharey=True,layout='constrained')
	ax=ax.ravel()
	# loop through sims
	for i in range(nsims):
		sim = getSimulation(home+sims[i])
		d0 = sdfread(0)
		times = read_pkl('times')
		xi1,xi2,_=getConcentrationRatios(d0)
		du0 = duarr[i,0,:]
		du1 = duarr[i,1,:]
		gyroratio = (xi1/xi2)*(marr[1]/marr[0])*(qarr[0]/qarr[1])**2
		du0_rat = du0/gyroratio
		du1_rat = du1/(1/gyroratio)
		# plot original traces
		ax[i].plot(times/tcnorm,du0,'r-')
		ax[i].plot(times/tcnorm,du1,'b-')
		# plot modified traces
		ax[i].plot(times/tcnorm,du0_rat,color='pink',linestyle='--')
		ax[i].plot(times/tcnorm,du1_rat,color='purple',linestyle='--')
		ax[i].annotate(labels[i],xy=(0.5,0.95),xycoords='axes fraction',ha='center',va='top')
	# legend
	legend = fig.legend([r'$\Delta u_D$',r'$\Delta u_{He3}$',r'$\Delta u_{D}^\prime$',r'$\Delta u_{He3}^\prime$'],\
						loc='upper center',ncol=4,bbox_to_anchor=(0.5,1.1),borderpad=0.1)
	# ax[0].legend(,loc='best')
	# limits and formatting
	ax[0].set_ylim(-0.1,5)
	ax[0].set_xlim(0,4.0)
	ax[0].locator_params(axis='x',nbins=4)
	fig.supxlabel(r'$t/\tau_{cp}$',**tnrfont)
	fig.supylabel(r'$\Delta u$'+' '+r'$[Jm^{-3}]$',**tnrfont)
	fig.savefig(home+'du_prime_all.pdf',bbox_inches='tight')
	pass

# plots the gyro-resonance theory vs emp. for values at multiple times (not through time, see PlotGyroResonance)
def gyro_time_compare(home,sims,identify_indmat,identify_markersmat,species=['Deuterons','He3','Protons'],mean_to=10,multipanel=True,\
					figname='multiple_times',plot_du=True,lims=(-1,22)):
	"""
		Plots the gyro-resonance for multiple times either in a multi panel or in a (default) singular plot. 
		Specify which markers to use and at which time indices. Can either plot du1/du2 ratio vs theory (1:1 line)
		or dE1/dE2 ratio against xi2 which tends to (m2/m1)*(q1/q2)**2 line.
		IN:
			home				: str, location of home dir of sims
			sims				: lst, list of sims in home directory
			identify_indmat 	: lst of int, integer indices of points to identify
			identify_markersmat : lst of str, types of markers used to identify each identified point
			species				: lst of str, name of species in order (maj1, maj2, min)
			mean_to				: int, mean to take energy to find diff
			multipanel			: bool, true for multiple panels for each sim in sims
			figname				: str, unique name of figure
			plot_du				: bool, True plots du ratio, False plots dE ratio vs xi2 (not time, see PlotGyroResonance)
			lims				: arr, the limits of the plot of gyro-resonance (limits symmetric in x and y)
		OUT:
			None
	"""
	getSimulation(home+sims[0])
	times = read_pkl('times')*getCyclotronFreq(sdfread(0),species[-1])/(2*const.PI) # normalised
	labels = ["{:.2f}".format(times[i]) for i in identify_ind[0]]
	m1 = getMass(species[0]) ; m2 = getMass(species[1]) 
	q1 = getChargeNum(species[0]) ; q2 = getChargeNum(species[1])
	if multipanel:
		fig,ax=plt.subplots(nrows=1+(len(sims)-1)//4,ncols=4,figsize=(12,6),sharex=True,sharey=True)
		fig.subplots_adjust(hspace=0.075,wspace=0.075)
		ax=ax.ravel()
		colors = ['k']*len(sims)
		alphas = [1]*len(identify_ind[0])
	else:
		fig,ax=plt.subplots(figsize=(6,4))
		if plot_du:
			ax.plot([0,10],[0,10],color='darkgrey',linestyle='--',linewidth=0.5,zorder=0)
		else:
			ax.axhline((m2/m1)*(q1/q2)**2,color='darkgrey',linestyle='--',zorder=0)
		colors = plt.cm.rainbow(np.linspace(0,1,len(sims)))
		alphas = [1]*len(identify_ind[0]) #np.linspace(0.5,1,len(identify_ind[0]))
	# loop through each sim (concentration)
	for i in range(len(sims)):
		print(i)
		getSimulation(home+sims[i])
		xi1,xi2,_ = getConcentrationRatios(sdfread(0))
		times=read_pkl('times')
		# load energy densities
		u1 = read_pkl(species[0]+'_KEdens')
		u2 = read_pkl(species[1]+'_KEdens')
		du1_du2 = (u1-np.mean(u1[:mean_to]))/(u2-np.mean(u2[:mean_to]))
		for j in range(len(identify_ind[i])): # time
			if multipanel:
				axj = ax[j]
				ax[i].annotate(r'$t/\tau_{cp}=$'+labels[i],xy=(0.05,0.85),xycoords='axes fraction',fontsize=16,fontname='Times New Roman')
			else:
				axj = ax
				if plot_du:
					axj.scatter([(xi1/xi2)*(m2/m1)*(q1/q2)**2],du1_du2[identify_ind[i][j]],marker=identify_markers[i][j],\
								zorder=1,facecolor=colors[i],edgecolor='none',s=30,alpha=alphas[j])
				else: # plot kinetic energy per-particle
					axj.scatter([xi2],(xi2/xi1)*du1_du2[identify_ind[i][j]],marker=identify_markers[i][j],\
								zorder=1,facecolor=colors[i],edgecolor='none',s=30,alpha=alphas[j])
		# plt.axhline((const.me_to_He3/const.me_to_mD)*(1/4),color='k')
		os.chdir('..')
	if multipanel: 
		xoff = 0.075 ; yoff = -0.01
		fig.supylabel(r'$\Delta u_D(t)/\Delta u_{He3}(t)$',x=xoff,**tnrfont)
		fig.supxlabel(r'$(\xi_{D}/\xi_{He3})(m_{He3}/m_D)(q_D/q_{He3})^2$',y=yoff,**tnrfont)
		figname = 'multipanel'
		ax.legend(labels,loc='best')
	else:
		if plot_du:
			# xoff = 0.05 ; yoff = -0.05
			ax.set_ylabel(r'$\Delta u_D(t)/\Delta u_{He3}(t)$',**tnrfont)
			ax.set_xlabel(r'$(\xi_{D}/\xi_{He3})(m_{He3}/m_D)(q_D/q_{He3})^2$',**tnrfont)
			figname = 'singlepanel_du'
			ax.set_xlim(lims) ; ax.set_ylim(lims)
		else:
			ax.set_ylabel(r'$\Delta E_D(t)/\Delta E_{He3}(t)$',**tnrfont)
			ax.set_xlabel(r'$\xi_{He3}$',**tnrfont)
			figname = 'singlepanel_dE'
		# axcopy = copy.copy(ax)
		# axinset = zoomed_inset_axes(ax, 0.25, loc=1)
		# axinset.add_artist(axcopy)
		# axinset.set_xlim(0.05,0.4)
		# axinset.set_ylim(-0.15,0.4)
	# plt.xlim(0,7)
	# plt.ylim(-1,22)
	# plt.show() ; sys.exit()
	fig.savefig('gyro_resonance_{}.png'.format(figname),bbox_inches='tight')
	return None

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
	
	# # D-He3
	# home = '/storage/space2/phrmsf/lowres_D_He3/'
	# sims = np.sort([i for i in os.listdir(home) if 'p_90' in i])[1:]
	# print(sims)
	# hlabels = np.array([int(i[2:4]) for i in sims])
	# # standard
	# xlabel=r'$(\xi_{D}/\xi_{He3})(m_{He3}/m_{D})(q_{D}/q_{He3})^2$'
	# ylabel=r'$[\Delta u_{D}/\Delta u_{He3}]_{max}$'
	# # # through-time
	# # ylabel=r'$(\Delta u_{D}/\Delta u_{He3})$' + ' ' + r'$(\xi_{He3}/\xi_D)$'
	# duarr = majIons_edens_ratio(home,sims,species=['Deuterons','He3'],norm_spec='Protons')
	# print(duarr,duarr.shape)
	# PlotGyroResonance(home,duarr,sims,species=['Deuterons','He3'],norm_spec='Protons',labels=hlabels,\
	# 				  ylabel=ylabel,xlabel=xlabel,lims=((0,8),(0,8)),through_time=False,identify=True,time_norm=r'$\tau_{cp}$')
	# # du_gyroratio(home,duarr,sims,species=['Deuterons','Helium3'],norm_spec='Protons',labels=hlabels)

	# p-B11
	home = '/storage/space2/phrmsf/'
	sims = ['p_B11']
	xlabel=r'$(\xi_{p}/\xi_{B11})(m_{B11}/m_{p})(q_{p}/q_{B11})^2$'
	ylabel=r'$[\Delta u_{p}/\Delta u_{B11}]_{max}$'
	duarr = majIons_edens_ratio(home,sims,species=['Protons','B11'],norm_spec='Alphas')
	PlotGyroResonance(home,duarr,sims,species=['Protons','B11'],norm_spec='Alphas',labels=[],xlabel=xlabel,ylabel=ylabel,\
						through_time=True,identify=False,time_norm=r'$\tau_{c\alpha}$')
	