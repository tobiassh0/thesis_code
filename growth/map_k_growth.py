
from func_load import *
#---------------------------------------------------------------------#
"""
def PITENSOR(file0,v0,kall,omegaall,theta=90):
	theta = theta*(const.PI/180) # radians
	PIxx = np.zeros(len(omegaall),complex) ; PIxy = np.zeros(len(omegaall),complex) ; PIyy = np.zeros(len(omegaall),complex)
	wci = getEffectiveCyclotronFreq(file0)

	for i in range(0, omegaall.shape[0]):
		l = round(omegaall[i]/wci) #l closest to the omega
		k = kall[i]
		kpara = kall[i]*np.cos(theta)
		kperp = kall[i]*np.sin(theta)
		
		za = kperp*v0/wci
		zarr = np.linspace(0,2*za,10000)
		dzarr = (zarr[-1]-zarr[0])/len(zarr)

		# try and except clause to remove l=0 division
		try:
			PIxx[i] = spec.jv(2*l,2*za)
		except:
			None
		try:
			PIxy[i] = (-1j*za/l)*(spec.jvp(2*l,2*za))
		except:
			None
		try:
			PIyy[i] = (1-(za/l)**2)*spec.jv(2*l,2*za)+(za/(2*l**2))*np.sum(dzarr*spec.jv(2*l,zarr))
		except:
			None

	return PIxx, PIxy, PIyy

def CHI0CALC(PIxx,PIxy,PIyy,n0,kall,omegaall,M,L,Z1,Z2,Z3,species=['Deuterons','Tritons','Alphas'],B=2.1):
	meff = (getMass(species[0])/Z1)*(1-M*Z2-L*Z3) + getMass(species[1])*M + getMass(species[2])*L
	print(PIxx.shape,PIxy.shape,PIyy.shape,meff.shape,omegaall.shape,kall.shape)
	Chi02 = PIxx + (const.qe*const.mu0*n0/(kall))*(((-2*1j)*omegaall/B)*PIxy + (1/meff)*PIyy)
	return np.sqrt(Chi02)


sim_loc = getSimulation('/storage/space2/phrmsf/traceT_D_99_T_01')
nval = 100
Mu = np.linspace(0,1,100)
Lambda = np.linspace(1e-5,1e-2,100)
Z1 = 1 ; Z2 = 1 ; Z3 = 2
CHI0 = np.zeros((len(Mu),len(Lambda),10000))
species = ['Deuterons','Tritons','Alphas']
B0 = 2.1
E0 = 3.5E6 * const.qe
v0 = np.sqrt(2*E0/getMass('Alphas'))
n0 = getMeanquantity(sdfread(0),'Derived_Number_Density_Electrons')

for m in range(len(Mu)):
	for l in range(len(Lambda)):
		meff = (getMass(species[0])/Z1)*(1-Mu[m]*Z2-Lambda[l]*Z3) + getMass(species[1])*Mu[m] + getMass(species[2])*Lambda[l]
		wci = const.qe*B0/meff
		omegas = wci*np.linspace(0,25,10*nval)
		_,k2,_ = coldplasmadispersion(sdfread(0),omegas)
		PIxx, PIxy, PIyy = PITENSOR(sdfread(0),v0,k2,omegas,theta=89)
		CHI0[m,l,:] = CHI0CALC(PIxx,PIxy,PIyy,n0,k2,omegas,Mu[m],Lambda[l],Z1=1,Z2=1,Z3=2,species=['Deuterons','Tritons','Alphas'],B=2.1)
dumpfiles(CHI0,'CHI0_Mu_Lambda')
plt.imshow(CHI0,**kwargs) ; plt.show()
"""

#---------------------------------------------------------------------#

# calculate the growth of k-modes for various time spacings
def growth_wavenumber(sim_loc,normspecies='Alphas',omega_min=0,omega_max=50,domega=0.25,dt_frac=0.5,\
					tstart_frac=0.5,tend_frac=2.0,theta=89.,color='k',plot=False):
	# define sim
	sim_loc = getSimulation(sim_loc)
	d0 = sdfread(0)
	
	# normalisation
	vA = getAlfvenVel(d0)	
	wnorm = getCyclotronFreq(d0,normspecies)
	klim = 0.5*2*const.PI/getdxyz(d0)
	times = read_pkl('times') ; tnorm = 2*const.PI/wnorm
	T = times[-1]/tnorm
	
	# cold plasma disp
	omegas = wnorm*np.arange(int(omega_min),int(omega_max),domega)
	species = getIonSpecies(d0)	
	_,k2,_ = coldplasmadispersion(d0, omegas, theta=theta)
	#k2 = omegas/(vA*np.sqrt(1+np.cos(86.3*const.PI/180.)**2)) # simpler dispersion
	thresh = k2 > 0
	k2 = k2[thresh] ; omegas = omegas[thresh]
	
	# FT1d & extract data at karr points
	FT1d = read_pkl('FT_1d_Magnetic_Field_Bz')
	harmonics = FT1d.shape[1]*(k2/klim)
	
	# loop through harmonics, then loop through time window to acquire multiple growth rates for a given k-mode
	[t0,t1] = [tstart_frac*tnorm,tend_frac*tnorm]
	t1_2 = (t1+t0)/2
	dt = dt_frac*tnorm
	nmin = 1
	nmax = (t1+t0)/(2*dt) ; narr = np.arange(nmin,nmax+1,1)
	growthRates = np.zeros((len(harmonics),len(narr)))	
	
	# loop through each harmonic index
	for i in range(len(harmonics)):
		# data = np.mean((FT1d[1:,round(harmonics[i])-1:round(harmonics[i])+1]),axis=1)
		print(round(harmonics[i]))
		data = (FT1d[1:,round(harmonics[i])])
		
		# extract peaks in data
		peaks = extractPeaks(data,Nperw=10)
		datapeak = data[peaks]
		time_dummy = times[peaks]
		# indexpeak = np.ones(len(datapeak))*i
		# axs.scatter(indexpeak,time_dummy/tnorm,datapeak)
		
		# fit data over multiple times
		for n in range(len(narr)):
			try:
				thresh = (time_dummy > t1_2-narr[n]*dt) & (time_dummy < t1_2+narr[n]*dt)
				timefit = time_dummy[thresh]
				datafit = datapeak[thresh]
				# assumed exponential growth (as expected from MCI)
				popt, pcov = curve_fit(lambda t,b,c: np.exp(b*t)+c,timefit,datafit,maxfev=5000) # exponential fitting
				b,c = popt
				fitted = np.exp(b*timefit) + c
				growthRates[i,n] = 2*const.PI*b # not normalised
			except:
				continue
	# === # 

	# mean growth rate for of all time ranges
	growthRatesMean = np.mean(growthRates,axis=1)
	growthRatesSTD = np.std(growthRates,axis=1)

	# plot
	if plot: # boolean
		thresh = growthRatesMean > 0
		growthRatesMean = growthRatesMean[thresh]
		omegas = omegas[thresh]
		plt.plot(omegas/wnorm,growthRatesMean/wnorm,'o-',color=color)
		plt.savefig('growth_kmodes.png')
	
	return omegas, growthRatesMean, growthRatesSTD

## multiple empirical growth rates (for multiple time-ranges) plotted against theory 
def multi_empirical_growths(home,sim_lst,labels,minions='Alphas',majions='Deuterons',maj2ions='Tritons',Emin=3.5e6,colors=[None],\
							theory_sim=None,theta_deg=89,times=[[0,0.5],[0.5,2.0]],wmin=0,wmax=25,domega=0.5,nval=int(1e2)):
	"""
		Function for an n-rows plot of the [0] theoretical growth rates [1] time range 1 [2] time range 2 ... [n] time range n
		Each time range in the times array is given as 2 values, times[n] = [tstart, tend]. Using these ranges, the growth rates
		are extracted empirically according to the cold plasma dispersion
		
		In:
			sim_lst 		: list of simulations to compare
			labels 		: (not used) can label each species 
			majions 		: majority ions
			minions 		: minority ions
			Emin 			: Minority species energy in eV
			theory_sim 	: simulation used to calculate theoretical growth rates
			theta_deg	: angle between B-field and sim domain (degrees)
			times 		: time array shape (n,2) which dictates time ranges to calculate empirical growth rates
			wmin & wmax : the minimum and maximum frequencies used for plotting & growth rates
			domega				: spacing in normalised freq to empirically calculate the growth rates
			nval			: the number of data points to calculate the theoretical growth rates
		Out:
			n-rows plot of the theory, and n time ranges		
	"""	
	## setup
	times = np.array(times) # convert to numpy array
	if len(colors) < len(sim_lst):
		colors = plt.cm.rainbow(np.linspace(0,1,len(sims)))
	shapes = ['o']*len(sim_lst)
	
	## plot setup
	fig,ax=plt.subplots(figsize=(8.5,4))
	ax_emp = ax # ax.twinx()
	
	## majority harmonics
	for i in np.arange(0,wmax+1,1):
		ax.axvline(i,color='darkgrey',linestyle='--')
	
	## theory
	_=getSimulation(theory_sim)
	d0 = sdfread(0)

	## setup	
	wcyca = getCyclotronFreq(d0,minions)
	wnorm = wcyca
	omegaall = wnorm*np.linspace(wmin,wmax,int(nval))
	vA = getAlfvenVel(d0)
	print(vA)

	## cold plasma dispersion
	_,kall,_ = coldplasmadispersion(d0,omegaall)
	v0 = np.sqrt(2*Emin*const.qe/getMass(minions))
	u = 0.98*vA #v0*np.sin(pitch_angle)
	vd = np.sqrt(v0**2 - u**2) #v0*np.cos(pitch_angle)
	vr = v0/1000
	
	## get theoretical growth rates & plot
	# posomega, posgamma = growth_rate_theory(minions, majions, theta, d0, u, vd, vr, kall, omegaall)
	posomega, posgamma = growth_rate_manual(minions=minions,majions=majions,maj2ions=maj2ions,wmax=wmax,\
											theta_deg=theta_deg,xi3=10**(-4),xi2=0.0,nval=int(nval))
	# get peak frequencies and plot as scatter
	peaks_th = extractPeaks(posgamma,Nperw=1)
	ax.scatter(posomega[peaks_th]/wcyca,posgamma[peaks_th]/wcyca,color='k',marker='x')
	ax.plot(posomega/wcyca,posgamma/wcyca,color='k')
	
	## empirical growths
	for ttimes in times:
		SimIndex=0
		for sim in sim_lst:

			## define sim
			sim_loc = getSimulation(home+sim)
			theta,_ = getMagneticAngle(sdfread(0))

			## get empirical growth rates for a range of k values and corresponding frequencies
			omegas, growthRatesMean, growthRatesSTD = growth_wavenumber(sim_loc,minions,wmin,wmax,domega,tstart_frac=ttimes[0],tend_frac=ttimes[1],theta=theta*180/const.PI)

			## thresh so if < 0 ignore data point
			thresh = growthRatesMean > 0
			growthRatesMean = growthRatesMean[thresh]
			omegas = omegas[thresh]

			## plot
			# ax_emp.plot(omegas/wnorm,growthRatesMean/wnorm,'-o',color=colors[SimIndex])
			ax_emp.scatter(omegas/wnorm,growthRatesMean/wnorm,color=colors[SimIndex],label=labels[SimIndex])
			SimIndex+=1
		ax_emp.set_xlim(0,wmax)
		

	## formatting
	# yscale
	ax.set_yscale('log')
	ax_emp.set_yscale('log')
	# ylabel
	# ax_emp.set_ylabel(r'$\gamma_s/$'+getOmegaLabel(minions),**tnrfont) # what colour to set this to?
	ax.set_ylabel(r'$\gamma/$'+getOmegaLabel(minions),**tnrfont)
	# xlabel
	ax.set_xlabel(r'$\omega/$'+getOmegaLabel(minions),**tnrfont)
	# # ylim
	ax.set_ylim(1e-2,10**(1)) # TODO; hardcoded
	# ax_emp.set_ylim(1e-6,1e1) # TODO; hardcoded
	# xlim
	ax.set_xlim(wmin,wmax)
	# legend
	ax_emp.legend(loc='best')

	# ## annotation label (alpha) 
	# ax_emp.annotate(str(ttimes[0])+r'$<t/\tau_{c\alpha}<$'+str(ttimes[1]),xy=(0.04,0.7125),xycoords='axes fraction',\
	# 					**tnrfont,ha='left',va='bottom',bbox=dict(boxstyle='square',pad=0.15,fc='w', ec='k', lw=1))
	
	fig.savefig(home+'referee_reports/Bz_kt_growth_rates.png',bbox_inches='tight')
	plt.show()
	
	# fig.savefig(home+'/Bz_kt_{}_{}_{}.png'.format(wmin,wmax,domega),bbox_inches='tight')
	# fig.savefig(home+'/Bz_kt_{}_{}_{}.eps'.format(wmin,wmax,domega),bbox_inches='tight')
	return None

#---------------------------------------------------------------------#

if __name__=='__main__':
	from func_load import * 
	# D-T
	home = '/storage/space2/phrmsf/traceT/'
	sims = ['traceT_D_100_T_00','traceT_D_99_T_01','traceT_D_89_T_11']
	tlabels = [r'$0\%$',r'$1\%$',r'$11\%$']
	colors = np.flip(['darkturquoise','r','g'])
	multi_empirical_growths(home,np.flip(sims),np.flip(tlabels),theory_sim=home+sims[0],\
							colors=colors,times=[[0.5,2.0]],nval=int(2e6),domega=0.25) # times=[[0.5,2.0]]
	# nval=int(2e6)
	
	# # D-He3
	# home = '/storage/space2/phrmsf/lowres_D_He3/'
	# wcp = const.qe*3.7/getMass('Protons')
	# # sims = np.sort([i for i in os.listdir(home) if 'p_90' in i])
	# # hlabels = np.array([int(i[2:4]) for i in sims])	
	# # multi_empirical_growths(home,sims,hlabels,'Deuterons','Protons',theory_sim=sims[0])
	# sims = ['0_00_p_90','0_45_p_90']
	# colors=['b','r']
	# # get empirical growth rates for a range of k values and corresponding frequencies
	# for i in range(len(sims)):
	# 	simloc=getSimulation(home+sims[i])
	# 	omegas, growthRatesMean, growthRatesSTD = growth_wavenumber(simloc,'Protons',0,20,domega=0.01,\
	# 												tstart_frac=0,tend_frac=5,theta=89*180/const.PI)
	# 	thresh = growthRatesMean > 0
	# 	plt.scatter(omegas[thresh]/wcp,growthRatesMean[thresh]/wcp,color=colors[i])
	# plt.show()
