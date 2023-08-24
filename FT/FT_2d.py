
from func_load import *


#sim_lst = ['traceT_0_00','traceT_D_99_T_01','traceT_D_89_T_11','traceT_0_50']
#sim_lst = ['D_He3_0_10_min_p_0_9']

sim_lst = ['ECRH/ECRH_JT60U_3']
home = '/storage/space2/phrmsf/'

for sim in (sim_lst):
	## setup
	fig, ax = plt.subplots(figsize=(8,4)) #(8,4)
	sim_loc = getSimulation(home+sim)
	d0 = sdfread(0)
	FT_2d = read_pkl('FT_2d_Magnetic_Field_Bz')
	times = read_pkl('times')
	dt = times[-1]/len(times)
	dx = getdxyz(d0)
	LDe = getDebyeLength(d0,'Electrons')
	## freq limits
	klim = 2*0.5*const.PI/dx
	wlim = 2*0.5*const.PI/dt
	vA = getAlfvenVel(d0)
	wcyc = getCyclotronFreq(d0,'Electrons')
	wcycp= getCyclotronFreq(d0,'Protons')
	wnorm = wcyc
	knorm = wcyc/vA
	klim_prime = klim/knorm
	wlim_prime = wlim/wnorm
	print(klim_prime,wlim_prime)
	(nw,nk) = FT_2d.shape
#	kmax = 40; wmax = 100
	kmax = klim_prime; wmax = wlim_prime
	dk = 2*knorm/(2*nk) ; dw = (2*wnorm/len(times))

#	_,_ = power(klim_prime=klim_prime,wlim_prime=wlim_prime,wmax=35,kmax=kmax,quantity='Magnetic_Field_DeltaBz',plot=False,read=False)

	# print dw, dk, klim and wlim in various units
	print('klim (vA/omega_p, va/omega_e) :: ',klim*vA/wcycp,klim/knorm) 
	print('wlim (omega_p, omega_e) :: ',wlim/wcycp,wlim/wnorm)
	print('dk (vA/omega_p, va/omega_e) :: ',dk*vA/wcycp,dk/knorm) 
	print('dw (omega_p, omega_e) :: ',dw/wcycp,dw/wnorm) 	

	# chopping and plotting
	FT_2d = FT_2d[:int(nw*wmax/wlim_prime),:int(nk*kmax/klim_prime)]
	ax.imshow(np.log10(FT_2d),interpolation='nearest',cmap='magma',origin='lower',aspect='auto',extent=[0,kmax,0,wmax],vmin=-4,vmax=6)
	ax.set_xlabel(r'$kv_A/\Omega_e$',fontsize=18) # v_A/\Omega_D
	ax.set_ylabel(r'$\omega/\Omega_e$',fontsize=18)

#	karea = 0.06
#	warea = 40
#	Nt = len(times)
#	nfft = Nt//5
#	noverlap = nfft//2
#	quantity = 'Magnetic_Field_Bz'
#	try: # dumped in one
#		fieldmatrix = read_pkl('fieldmatrix_'+quantity)
#	except: # batch dumped
#		fieldmatrix = load_batch_fieldmatrix(list_sdf(sim_loc),quantity,para=False)
##	T = times[-1]
##	L = getGridlen(d0)
##	_, _ = getBicoh(karea,warea,fieldmatrix,dt,T,L,wcyc,1/LDe,nfft=nfft,noverlap=noverlap,window=True,bispectrum=True)		

	## ions and characteristic freq
	maj_species, maj2_species, min_species = getIonSpecies(d0)
	wce  = getCyclotronFreq(d0,'Electrons',Z=1)
	wpe	 = getPlasmaFreq(d0,species='Electrons')
	wpi  = getPlasmaFreq(d0,species='Protons')

	file0 = d0
#	nmaj1 = getMeanquantity(file0,'Derived_Number_Density_'+maj_species)
#	n_e = getMeanquantity(file0,'Derived_Number_Density_Electrons')
#	nmaj2 = 0
#	nmin = 0
#	if maj2_species != '':
#		nmaj2 = getMeanquantity(file0,'Derived_Number_Density_'+maj2_species)
#	if min_species != '': # check if thermal
#		nmin = getMeanquantity(file0, 'Derived_Number_Density_'+min_species)
#	Kappa = nmaj1/n_e ; Mu = nmaj2/n_e ; Alpha = nmaj2/nmaj1 ; Lambda = nmin/n_e ; Xi = nmin/nmaj1
#	print(Kappa,Mu)
#	m_eff = Kappa*getMass(maj_species)+Mu*getMass(maj2_species)
#	w_PI2 = (n_e*const.qe**2)/(m_eff*const.e0)
#	w_CI2 = (const.qe*getMeanField3D(d0,'Magnetic_Field_B')/m_eff)**2
#	W_LH = (((w_PI2)+(w_CI2))/(1+((wpe**2)/(wce**2))))**.5	
#	ax.axhline(W_LH/wnorm,color='white',linestyle='-')
#	print(W_LH/wnorm)
	
	## debugging
	print(wpi)
	print(wpe)
	print(wcyc)
	print(wce)
	print(maj_species,maj2_species,min_species)

# uncomment for cold-plasma dispersion overplotted	
	## cold plasma disp & Wave modes
	omegas = wnorm*np.linspace(0,wmax,10000)
	k1,k2,k3=coldplasmadispersion(d0,maj_species,maj2_species,omegas=omegas)
	ax = ColdWaveModes(ax,[wpi,wpe,wcyc,wce],vA,wnorm,LH=True,W2=True)
	k1,k2,k3 = k1/knorm, k2/knorm, k3/knorm 
	thresh = k2 > 0 
	ax.plot(k2[thresh],omegas[thresh]/wnorm,color='k',linestyle='-',alpha=0.75)
	ax.set_ylim(0,wmax)
	ax.set_xlim(0,kmax)
	plt.show()
#	fig.savefig(home+sim+'_FT_2d_Bz_v2.png',bbox_inches='tight')
