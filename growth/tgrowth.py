
def growth_rate_manual(minspec='Alphas',majions='Deuterons',maj2ions='Tritons',wmax=25,theta_deg=89,xi3=10**(-4),xi2=0.0,\
					B0=2.1,ne=1e19,pitchangle=0.69115,_Emin=3.5,vspread=1/100,nval=int(1e5),plot=False):
	"""
	Plots the linear growth rate, from eqs. (29) and (8) in Refs. (R. Dendy et al., 1994) and (K. G. McClements et al., 1995) respectively. 
	Doesn't load a sim value, instead uses user-defined inputs for the initial conditions and plasma freqs.
		In:
			minspec 		: minority species
			majions 		: 1st bulk ion species
			maj2ions		: 2nd bulk ion species
			wmax 			: maximum frequency, as normalised to the minspec cyclotron freq 
			theta_deg 	: angle between magnetic field and x-hat (in degrees)
			xi3			: minority species concentration
			xi2			: 2nd bulk ion species concentration
			B0				: magnetic field strength (Tesla)
			ne				: electron number density
			pitchangle	: the cosine of the pitch angle, default of 0.22*PI radians (as opposed to -0.646)
			_Emin			: the energy of the minority species (in units of MeV)
			vspread		: the thermal spread (vr/v3) of the parallel component 
			nval			: number of valuation points to calculate (aim for dw < 0.01 Omega_c,minspec)
			plot			: boolean to determine whether to plot data
		Out:
			omega 	: frequencies calculated for positive growth rates 
			gamma		: positive growth rates (gamma<0 replaced with 0)
	"""
	# initial conditions
	# B0 = 2.1
	# ne = 1e19
	theta = theta_deg*(const.PI/180) # rad
	del theta_deg
	
	# species masses, charges and densities
	Z1 = getChargeNum(majions) ; Z2 = getChargeNum(maj2ions); Z3 = getChargeNum(minions)
	# majions = 'Deuterons'
	# maj2ions= 'Tritons' 
	# minions = 'Alphas'
	m1 = getMass(majions)
	m2 = getMass(maj2ions)
	m3 = getMass(minions)
	n1 = ne * (1/Z1)*(1-Z2*xi2-Z3*xi3)
	n2 = xi2*ne
	n3 = xi3*ne
	
	# checking quasi-neutrality
	if ne != (n1*Z1 + n2*Z2 + n3*Z3): # @assert quasi-neutrality
		print('Quasi-neutrality not maintained')
		raise SystemExit
	else:
		print('Quasi-neutrality maintained')
	
	# Alfven velocity
	rho0 = ne*const.me + n1*m1 + n2*m2 + n3*m3
	vA = B0/np.sqrt(const.mu0*rho0)
	
	# species frequencies
	wca = Z3*const.qe*B0/m3
	wci = Z1*const.qe*B0/m1
	wce = const.qe*B0/const.me
	wpa = np.sqrt((n3*(Z3*const.qe)**2)/(m3*const.e0))
	wpi = np.sqrt((n1*(Z1*const.qe)**2)/(m1*const.e0))
	wpe = np.sqrt((ne*(const.qe)**2)/(const.me*const.e0))

	## frequencies can be defined in two ways
	# (1) cold plasma dispersion
	omegas = wca*np.linspace(0,wmax,2*nval)	
	k1, k2, k3 = coldplasmadispersion_analytical(omegas,wpf=[wpe,wpa,wpi],wcf=[wce,wca,wci],theta=theta)

	# (2) freq of EM wave, eq. () from ()
	# k2 = (wca/vA)*np.linspace(0,20,nval)
	# kpara = k2 * np.cos(theta)
	# kperp = k2 * np.sin(theta)
	# omegas = (0.5*vA**2)*(k2**2 + kpara**2 + (k2*kpara*vA/wci)**2 + ((k2**2 + kpara**2 + (k2*kpara*vA/wci)**2)**2 - (2*k2*kpara)**2)**0.5)

	# velocity and energies
	Emin = (_Emin * 10**6)*const.qe # minority energy
	v3 = np.sqrt(2*Emin/m3) # minority velocity
	# u_vA = 0.98 # uperp/vA
	# u = u_vA * vA # perp drift
	# u = v3 * np.cos(-0.646)
	u = v3 * np.cos(pitchangle) # as per traceT sims
	vd = np.sqrt(v3**2 - u**2) # para drift
	vr = v3*vspread # para spread

	# wavenumbers
	kpara = k2 * np.cos(theta)
	kperp = k2 * np.sin(theta)
	Npara = kpara*vA/omegas
	Nperp = kperp*vA/omegas
#	larr = (omegas/wca)//1 # find nearest harmonic as array

	# growth rates
	gamma = np.zeros(len(omegas))
	for i in range(len(omegas)):
		l = np.int32(omegas[i]/wca) # find nearest harmonic
		eetal = (omegas[i]-kpara[i]*vd-l*wca)/(kpara[i]*vr)
		za = kperp[i]*u/wca
		Jl = spec.jv(l,za) # (order, argument)
		Jlprime = spec.jvp(l,za)
		JlJlprime = Jl*Jlprime
#		plt.scatter(za,Jl)

		# Ml term #
		ml1 = 2*l*(omegas[i]/wci)*((Jlprime**2) + (1/(za**2))*(l**2 - za**2)*(Jl**2))
		ml2 = -2*((omegas[i]**2 - wci**2)/(wci**2))*(JlJlprime/za)*((l**2)*(Nperp[i]**2)+(2*l**2 - za**2)*(Npara[i]**2))
		ml3 = 2*(JlJlprime/za)*(za**2 - 2*l**2)
		Ml = ml1 + ml2 + ml3
		
		# Nl terms #
		nl1 = -2*l*(omegas[i]/wci)*(JlJlprime/za) + ((omegas[i]**2 - wci**2)/((wci**2)))*(Npara[i]**2 * (((l**2 * Jl**2)/(za**2)) + Jlprime**2) + Nperp[i]**2 * ((l**2 * Jl**2)/(za**2)))
		nl2 = ((l**2 * Jl**2)/(za**2)) + Jlprime**2
		Nl = nl1 + nl2
		
		# gamma #
		g1 = ((wpa**2)/(wpi**2))*((wci**4)/((wci+(omegas[i]-wci)*Npara[i]**2)*(wci-(omegas[i]+wci)*Npara[i]**2)))
		g2 = ((l*wca/(kpara[i]*vr))*Ml - 2*((u**2)/(vr**2))*eetal*Nl)
		g3 = (np.sqrt(const.PI)/(2*omegas[i])) * np.exp(-1*(eetal**2))
		gamma[i] = g1 * g2 * g3
	
	if plot:		
		plt.plot(omegas/wca,np.log10(gamma/wca))		
		plt.ylabel('log10(gamma/wca)')
		plt.xlabel('w/wca')
		plt.ylim(-6,4)
		os.chdir('/storage/space2/phrmsf/traceT/growth_rates')
		plt.savefig('{}_{}_{}_{}_growthrates.png'.format(B0,ne,xi2,xi3))
		plt.clf()

	# positive growth rates
	gamma[gamma<0] = 0
	return omegas, gamma



if __name__=='__main__':
	from func_load import *
	for txi2 in [0.01,0.05,0.11,0.5]:
		omega, gamma = growth_rate_manual(xi2=txi2)
	sys.exit()

	# load example sim (densities etc.)
	sim_loc = getSimulation('/storage/space2/phrmsf/traceT/old/traceT_highres_0_01')
	ind = list_sdf(sim_loc)
	d0 = sdfread(0)
	
	# number of compute points
	nval = int(1e6)
	
	# species
	minions = 'Alphas'
	majions = 'Deuterons'
	elec = 'Electrons'
	
	# parameters
	theta = 89. # deg
	wc_maj = getCyclotronFreq(d0,majions)
	vA = getAlfvenVel(d0)
	E0 = 3.5E6 * const.eV_to_J
	malpha = const.me*const.me_to_malpha
	
	# birth and spread velocities
	v0 = np.sqrt(2*E0/malpha)
	u = np.cos(0.22*const.PI)*v0
	vd = np.sin(0.22*const.PI)*v0
	vr = 0.001*v0#*(1/(np.sqrt(2)))
	print('vd/v0 :: {},\nvr/v0 :: {}\nu/v_A :: {}'.format(vd/v0,vr/v0,u/vA))
	
	# setup freq range
	omegas = wc_maj*np.linspace(0,25,nval)
