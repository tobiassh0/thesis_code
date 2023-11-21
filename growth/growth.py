
def TGROWTH():
	sim_loc = getSimulation('/storage/space2/phrmsf/traceT/old/traceT_highres_0_01')
	ind = list_sdf(sim_loc)
	nval = int(1e6)
	
	sim_loc = getSimulation('/storage/space2/phrmsf/traceT/old/traceT_highres_0_01')
	theta = 89.
	minions = 'Alphas'
	majions = 'Deuterons'
	wcyca = getCyclotronFreq(sdfread(0),minions)
	omegaall = wcyca*np.linspace(0,25,nval)
	vA = getAlfvenVel(sdfread(0))
	print(vA)
	_,kall,_ = coldplasmadispersion(sdfread(0),omegaall)
	emin = 3.5E6
	v0 = np.sqrt(2*emin*const.qe/getMass(minions))
	val = 0.001
	u = 0.98*vA #v0*np.sin(pitch_angle)
	vd = np.sqrt(v0**2 - u**2) #v0*np.sin(pitch_angle)
	vr = v0/1000
	w2, gamma2 = growth_rate_man(minions, majions, theta, sdfread(0), u, vd, vr, kall, omegaall)
	plt.plot(w2/wcyca,gamma2/wcyca) ; plt.show()
	sys.exit()

	minions = 'Alphas'
	majions = 'Deuterons'
	elec = 'Electrons'
	theta,_ = getMagneticAngle(sdfread(0))
	wc_min = getCyclotronFreq(sdfread(0),minions)
	vA = getAlfvenVel(sdfread(0))
	E0 = 3.5E6 * const.eV_to_J
	malpha = const.me*const.me_to_malpha
	v0 = np.sqrt(2*E0/malpha)
	u = 0.98*vA # np.cos(0.22*const.PI)*v0
	vd = np.sqrt(v0**2 - u**2) # np.sin(0.22*const.PI)*v0 #0.06*v0#0.001*v0#*(1/(np.sqrt(2)))
	vr = v0/1000
	print('vr/v0 :: {}\nu/v_A :: {}'.format(vr/v0,u/vA))
	omegas = wc_min*np.linspace(0,25,nval)
	thetas = [89]
	
	fig,ax = plt.subplots(figsize=(7,3))
	_,k2,_ = coldplasmadispersion(sdfread(0),omegas)
	posomega, posgamma = growth_rate_man(minions, majions, thetas[0], sdfread(0), u, vd, vr, k2, omegas)	
	#posgamma = np.array(posgamma) #; posgamma[np.isnan(posgamma)] = 0 
	wnorm = wc_min
	ax.plot(posomega/wnorm, posgamma/wnorm)
	ax.annotate(r'$\theta=$'+str(np.around(thetas[0],1))+r'$^\circ$',xy=(0.1,0.8),xycoords='axes fraction',fontsize=18)
	ax.locator_params(axis='y',nbins=4)
	ax.locator_params(axis='x',nbins=10)
	ax.set_ylabel(r'$\gamma/\Omega_\alpha$',fontsize=18)	
	ax.set_xlabel(r'$\omega/\Omega_\alpha$',fontsize=18)
	plt.show()


if __name__=='__main__':
	from func_load import *
	TGROWTH()
	nval = int(1e6)
	theta = 89.
	n0 = 1e19
	xib = 1e-3
	nb = n0*xib	
	ni = n0 - 2*nb # quasi-neutrality
	if ni + 2*nb != n0:
		print('# ERROR # : quasi-neutrality not asserted')
		raise SystemExit
	B0 = 2.1
	rho0 = nb*getMass('Alphas') + ni*getMass('Deuterons')
	vA = B0/np.sqrt(const.mu0*rho0)
	emin = 3.5e6 * const.qe
	v0 = np.sqrt(2*emin/getMass('Alphas'))
	u = v0*np.cos(-0.646)
	print(u/v0,u/vA)
	wce = const.qe*B0/const.me
	wcycb = 2*const.qe*B0/getMass('Alphas')
	wcyci = const.qe*B0/getMass('Deuterons')
	wpe = np.sqrt((n0*(const.qe)**2)/(const.me*const.e0))
	wpb = np.sqrt((nb*(2*const.qe)**2)/(getMass('Alphas')*const.e0))
	wpi = np.sqrt((ni*(const.qe)**2)/(getMass('Deuterons')*const.e0))
	wcyc = [wcycb,wcyci]
	wp = [wpb,wpi]
	omegaall = wpb*np.linspace(0,25,nval)
	_,kall,_=coldplasmadispersion_analytical(omegaall,wpf=[wpe,wpb,wpi],wcf=[wce,wcycb,wcyci],theta=theta)
	plt.plot(kall,omegaall/wcycb) ; plt.show()
	wprime, gprime = growth_rates_analytical_all(vA,theta,v0,u,kall,omegaall,val=1/1000,wcyc=wcyc,wp=wp)
	plt.plot(wprime,gprime,color='k') ; plt.show()