
def growth_aneutronic(home,sims,labels):
	ls = len(sims)
	fig,axs=plt.subplots(nrows=3,ncols=ls//3,figsize=(8,5))
	# extract basic parameters (how to include second species?)
	for i in range(ls):
		simloc = getSimulation(home+sims[i])
		d0 = sdfread(0)
		species = getIonSpecies(d0)
		majions, maj2ions, minions = species # split species 
		theta_deg,_ = getMagneticAngle(d0)
		xi1,xi2,xi3 = getConcentrationRatios(d0)
		vA = getAlfvenVel(d0)
		print("{:.2f}".format(xi2),vA/const.c)

#		print(((100*xi2)//1)/100) # round to nearest 100th
#		posomega, posgamma = growth_rate_manual(minions=minions,majions=majions,maj2ions=maj2ions,wmax=wmax,theta_deg=theta_deg,xi3=10**(-4),xi2=)
#		axs[i].plot(posomega,posgamma)

	return None

if __name__=='__main__':
	from func_load import *
	import map_k_growth as gk

	# D-He3
	home = '/storage/space2/phrmsf/lowres_D_He3/'
	fig,ax=plt.subplots(figsize=(8,8),nrows=2,sharex=True,sharey=True)
	wcp = const.qe*3.7/getMass('Protons')
	times = read_pkl(home+'0_00_p_90/times')

	# === empirical === # 
	# sims = np.sort([i for i in os.listdir(home) if 'p_90' in i])
	# hlabels = np.array([int(i[2:4]) for i in sims])	
	# multi_empirical_growths(home,sims,hlabels,'Deuterons','Protons',theory_sim=sims[0])
	sims = ['0_00_p_90','0_45_p_90']
	colors=['b','r']
	dw = (2*const.PI/times[-1])/wcp
	# get empirical growth rates for a range of k values and corresponding frequencies
	for i in range(len(sims)):
		simloc=getSimulation(home+sims[i])
		for j in range(0,21): # plot integer harmonics
			ax[i].axvline(j,color='darkgrey',linestyle='--')
		omegas, growthRatesMean, growthRatesSTD = gk.growth_wavenumber(simloc,'Protons',0,20,domega=0.001,\
													tstart_frac=0,tend_frac=5,theta=89*180/const.PI)
		thresh = growthRatesMean > 0
		ax[i].scatter(omegas[thresh]/wcp,growthRatesMean[thresh]/wcp,color=colors[i],zorder=1)
	
	# === theoretical === # 
	ximin = 1e-3
	u0 = np.sqrt(2*14.68*1e6*const.qe/getMass('Protons'))
	uperp_vA = 0.9

	# all D
	rhoM = 5e19*((1-ximin)*getMass('Deuterons') + ximin*getMass('Protons'))
	vA = 3.7/np.sqrt(const.mu0*rhoM)
	upara = np.sqrt(u0**2 - (uperp_vA*vA)**2)
	pitch = np.arctan(uperp_vA*vA/upara) # radians
	omegas, gamma = growth_rate_manual(minions='Protons',majions='Deuterons',maj2ions='He3',wmin=1,wmax=20,theta_deg=89,xi3=ximin,\
						xi2=0.0,B0=3.7,ne=5e19,pitchangle=pitch,_Emin=14.68,vspread=1/100,nval=int(1e5),plot=False)
	ax[0].plot(omegas/wcp,gamma/wcp,color='b',alpha=0.75,label=r'$D$',zorder=2)

	# all He3
	xi2 = 0.5*(1-ximin) # maximal xi2 
	rhoM = 5e19*(xi2*getMass('He3') + ximin*getMass('Protons'))
	vA = 3.7/np.sqrt(const.mu0*rhoM)
	upara = np.sqrt(u0**2 - (uperp_vA*vA)**2)
	pitch = np.arctan(uperp_vA*vA/upara) # radians
	omegas, gamma = growth_rate_manual(minions='Protons',majions='He3',maj2ions='Deuterons',wmin=1,wmax=20,theta_deg=89,xi3=ximin,\
						xi2=0.0,B0=3.7,ne=5e19,pitchangle=pitch,_Emin=14.68,vspread=1/100,nval=int(1e5),plot=False)
	ax[1].plot(omegas/wcp,gamma/wcp,color='r',alpha=0.75,label=r'$^3He$',zorder=2)
	
	# formatting
	fig.supylabel(r'$\gamma/\Omega_p$',**tnrfont)
	ax[-1].set_xlabel(r'$\omega/\Omega_p$',**tnrfont)
	ax[0].set_title('D')
	ax[-1].set_title('He3')
	# ax[-1].legend(loc='best')
	fig.savefig(home+'growth_rate_theory_and_empirical.png',bbox_inches='tight')
	ax[0].set_xlim(0,20) ; ax[0].set_ylim(0,0.5)
	plt.show()
