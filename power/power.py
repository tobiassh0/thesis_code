
from func_load import *


def power(klim_prime,wlim_prime,wmax,kmax,norm_omega=r'$\Omega_D$',quant='Magnetic_Field_Bz',plot=False,read=True,outp=True):

	if read:
		try:
			log10_power = read_pkl('log10_power')
			omegas = read_pkl('omegas_power')
			read = False; calc = False
		except:
			print('Can\'t read power and omegas...')
			read = False ; calc = True
	else:
		calc = True
	if calc:
		try:
			FT_2d = read_pkl('FT_2d_'+quant)
		except:
			print('## ERROR ## :: Cannot load FT_2d')
			raise SystemExit
		print('Calculating power...')
		log10_power,omegas=powerspectrum(FT_2d,wlim_prime,klim_prime,0,wmax,0,kmax)
		dumpfiles(log10_power,'log10_power')
		dumpfiles(omegas,'omegas_power')
	
	if plot:
		print('Plotting Power...')
		width,height = 8,5
		fig,ax=plt.subplots(figsize=(width,height))

		ax.plot(omegas,10**log10_power)
		for i in range(1,int(wmax)+1):
			ax.axvline(i,color='k',alpha=0.5,linestyle=':')
		ax.set_xlabel(r'$\omega$'+'/'+norm_omega,fontsize=18)
		ax.set_xlim(0,wmax)
		ax.set_ylabel('Power',fontsize=18) ## unitless
		## find min, max
		paxmin, paxmax = min(10**log10_power[10:]), max(10**log10_power[10:])
		ax.set_ylim(paxmin/10,10*paxmax) ## power of 10 higher and lower than min / max
		ax.set_yscale('log')
		plotting(fig,ax,'power')

	if outp: 
		return omegas, log10_power
	else: 
		del omegas; log10_power ; return None


#sim_lst = ['traceT_0_00','traceT_0_01','traceT_0_11']
#wmax = 30
#kmax = 100
#for sim in sim_lst:
#	sim_loc = getSimulation('/storage/space2/phrmsf/'+sim)
#	d0 = sdfread(0)
#	nx = len(getQuantity1d(d0,'Electric_Field_Ex'))
#	L = getGridlen(d0)
#	times = read_pkl('times')
#	nt = len(times)
#	vA = getAlfvenVel(d0)
#	wcyc = getCyclotronFreq(d0,'Deuterons')
#	FT_2d = read_pkl('FT_2d_Magnetic_Field_Bz')
#	klim = 0.5*2*nx*const.PI/L
#	klim_prime = klim*vA/wcyc
#	wlim = 0.5*2*nt*const.PI/times[-1]
#	wlim_prime = wlim/wcyc
#	maj_species = 'Deuterons'
#	omegas, log10_power = power(klim_prime=klim_prime,wlim_prime=wlim_prime,wmax=wmax,kmax=kmax,norm_omega=getOmegaLabel(maj_species),plot=False,read=False)
#	plt.plot(omegas,log10_power)
#for i in range(0,wmax):
#	plt.axvline(i,color='k',linestyle='--',alpha=0.3)
#	plt.axvline(i*(const.me_to_mT/const.me_to_mD),color='k',linestyle='-.',alpha=0.3)
#plt.legend(loc='best',labels=sim_lst)
#plt.show()

