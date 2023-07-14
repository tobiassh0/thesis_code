from list_new import *
import my_constants as const

kwargs={'interpolation':'nearest','origin':'lower','aspect':'auto'}


def coldplasma_eff(wci,wpi,wce,wpe,omegas,theta=None):
	theta = theta*const.PI/180
	sin = np.sin(theta) ; cos = np.cos(theta)
	
	wpf = [wpe, wpi]
	wcf = [wce , wci]

	l = len(omegas)
	R = np.ones(l) ; P = np.ones(l) ; L = np.ones(l) ; S = np.zeros(l) ; D = np.zeros(l) 
	B = np.zeros(l) ; F = np.zeros(l) ; A = np.zeros(l) ; C = np.zeros(l) 
	
#	for i in range(len(wpf)):	
#		R = R - ((wpf[i]**2)/(omegas*(omegas + wcf[i])))
#		L = L - ((wpf[i]**2)/(omegas*(omegas - wcf[i])))
#		P = P -  ((wpf[i]**2)/omegas**2)

	R = R - ((wpf[0]**2)/(omegas*(omegas + wcf[0]))) - ((wpf[1]**2)/(omegas*(omegas + wcf[1])))
	L = L - ((wpf[0]**2)/(omegas*(omegas - wcf[0]))) - ((wpf[1]**2)/(omegas*(omegas - wcf[1])))
	P = P -  ((wpf[0]**2)/omegas**2) - ((wpf[1]**2)/omegas**2)

	S = 0.5*(R+L) ; D = 0.5*(R-L)
	C = P*R*L
	B = R*L*(sin**2) + P*S*(1.0 +cos**2)
	F = (((R*L - P*S)**2)*(sin**4) + 4.0*(P**2)*(D**2)*(cos**2))**0.5
	A = S*(sin**2) + P*(cos**2)
	n1 = np.zeros(l, dtype=complex) ; n2 = np.zeros(l, dtype=complex) ; n3 = np.zeros(l, dtype=complex) ; n4 = np.zeros(l, dtype=complex) 
	n3 = np.lib.scimath.sqrt((R*L)/S)
	n1 =  np.lib.scimath.sqrt((B+F)/(2.0*A))
	n2 = np.lib.scimath.sqrt((B-F)/(2.0*A))
	n3 = -np.lib.scimath.sqrt((B+F)/(2.0*A))
	n4 = -np.lib.scimath.sqrt((B-F)/(2.0*A))
	del R, P, L, S, D, B, F, A
#	return n1*omegas/const.c, n2*omegas/const.c, n3*omegas/const.c, n4*omegas/const.c
	return (np.real(n1)*omegas)/const.c , (np.real(n2)*omegas)/const.c , np.real((n3*omegas)/const.c) #, (n4*omegas)/c, omegas


## Find the growth rates of the MCI in its linear phase based off of drift and spread velocities
def growth_rate_Eff(minions, theta, wcmin, wci, wpmin, wpi, va, u, vd, vr, kall, omegaall):

	wcycb = wcmin
	wpb = wpmin
	wcyci = wci
	wpi = wpi
	va = va

	theta = theta*(const.PI/180.0) #radians
	gammas = np.zeros(omegaall.shape[0]) #growth rates

	for i in range(0, omegaall.shape[0]):
		l = round(omegaall[i]/wcycb) #l closest to the omega
		k = kall[i]
		kpara = kall[i]*np.cos(theta)
		kperp = kall[i]*np.sin(theta)

		Npara = (kpara*va)/omegaall[i]
		Nperp = (kperp*va)/omegaall[i]

		eetal = (omegaall[i] - kpara*vd - l*wcycb)/(kpara*vr) # vd=0
		za = kperp*u/wcycb
		############## M_l ###############################################################
		mlterm1 = 2.0*l*(omegaall[i]/wcyci)*((spec.jvp(l,za)**2) + ((1.0/za**2)*(l**2 - za**2)*spec.jv(l,za)**2))
		mlterm2 = -2.0*((omegaall[i]**2 - wcyci**2)/wcyci**2)*((spec.jv(l,za)*spec.jvp(l,za))/za)*((l**2)*Nperp**2 - (za**2 - 2.0*(l**2))*(Npara**2))
		mlterm3 = (2.0*spec.jv(l,za)*spec.jvp(l,za)/za)*(za**2 - 2.0*(l**2))
		ml = mlterm1 + mlterm2 + mlterm3
		##################################################################################
						#
		############# N_l #################################################################
		nlterm1 = -2.0*l*(omegaall[i]/wcyci)*(spec.jv(l,za)*spec.jvp(l,za)/za)
		nlterm2pre = (omegaall[i]**2 - wcyci**2)/wcyci**2
		nlterm2 = (Npara**2)*((l*spec.jv(l,za)/za)**2 + spec.jvp(l,za)**2) + (Nperp*l*spec.jv(l,za)/za)**2
		nlterm3 = (l*spec.jv(l,za)/za)**2 + spec.jvp(l,za)**2
		nl = nlterm1 + nlterm2pre*nlterm2 + nlterm3
		#################################################################################		
						#
		########## Gamma ################################################################
		pre = (((wpb*(wcyci**2))/wpi)**2)*((const.PI**0.5)/(2.0*omegaall[i]))*np.exp(-(eetal**2))
		term1 = 1.0/((wcyci + (omegaall[i] - wcyci)*(Npara**2))*(wcyci - (omegaall[i] + wcyci)*(Npara**2)))
		term2 = ((l*wcycb*ml)/(kpara*vr)) - (2.0*eetal*nl)*((u/vr)**2)
		gammas[i] = pre*term1*term2
		#print gammas[i]
	
	####### Want gamma > 1 only ###########
	posomega = [] ; posgamma = []	

	for i in range(0,gammas.shape[0]):
		posomega.append(omegaall[i])
		if (gammas[i] >= 0):
			posgamma.append(gammas[i])
		else:
			posgamma.append(0)

	return np.array(posomega,dtype='float'), np.array(posgamma,dtype='float')
	




fig,axs=plt.subplots(nrows=1,ncols=1,figsize=(8,6))
B0 = 2.1 # background B-field
theta = 89. ## background B-field angle
maj1spec = 'Deuterons' ; maj2spec = 'Tritons' ; minspec = 'Alphas' ## species names
Z0 = getChargeNum('Electrons') ; Z1 = getChargeNum(maj1spec) ; Z2 = getChargeNum(maj2spec) ; Z3 = getChargeNum(minspec) ## real charge numbers
m0 = getMass('Electrons'); m1 = getMass(maj1spec) ; m2 = getMass(maj2spec) ; m3 = getMass(minspec) ## real masses
n0 = 1e19 ## electron number density
## ring beam parameters ##
E0 = 3.5E6 * const.eV_to_J
malpha = const.me*const.me_to_malpha
v0 = np.sqrt(2*E0/malpha)
u = np.cos(0.22*const.PI)*v0
vd = np.sin(0.22*const.PI)*v0
vr = 0.001*v0
wc1 = Z1*const.qe*B0/m1 ; wc3 = Z3*const.qe*B0/m3
##	##
xi3 = 1e-3 ## min particle concentration
## Parameters : 
#### zoom:  0 < omega/wc1 < 0.5
####		0.7 < xi2 < 1.0
#### full:	0 < omega/wc1 < 35
####		0 < xi2 < 1
xi2 = np.linspace(0.7,1,4000) # np.arange(0,1,0.0005)
omegaall = wc1*np.linspace(0,0.5,250000) # wc1*np.arange(0,35,0.0002) ## range of frequencies to calculate over
gamma = np.zeros((len(xi2),len(omegaall)))

## parallel
def PARAGROWTH(index):
	print(100*index/len(xi2))
	_,k2eff,_=coldplasma_eff(WCeff[index],WPeff[index],wce,wpe,omegaall,theta)
	posomegaEff,posgammaEff = growth_rate_Eff(minspec, theta, wc3, WCeff[index], wp3, WPeff[index], va[index], u, vd, vr, k2eff, omegaall)
	return posgammaEff

## arrays of values
xi1 = (1/Z1)*(1-Z2*xi2-Z3*xi3)
meff = (xi1*m1+xi2*m2+xi3*m3)
WCeff = const.qe*B0/meff
WPeff = np.sqrt((n0*(const.qe)**2)/(meff*const.e0))
wce = const.qe*B0/m0
wpe = ((n0*(const.qe**2))/(m0*const.e0))**0.5
wp3 = np.sqrt((n0*xi3*(Z3*const.qe)**2)/(m3*const.e0))
va = B0/np.sqrt((meff*n0)*const.mu0)

try:
	gamma = read_pkl('GAMMA_HEAT_zoom') # read_pkl('GAMMA_HEAT') # read_pkl('GAMMA_HEAT_zoom')
except:
	pool=mp.Pool(mp.cpu_count()-2)
	gamma=np.vstack(np.array(pool.map_async(PARAGROWTH,np.arange(0,len(xi2),1)).get(99999)))
	pool.close()
	dumpfiles(gamma,'GAMMA_HEAT_zoom') # read_pkl('GAMMA_HEAT') # read_pkl('GAMMA_HEAT_zoom')

	
im = axs.imshow(np.log10(gamma+1e-50).T,**kwargs,cmap='jet',extent=[xi2[0],xi2[-1],omegaall[0]/wc1,omegaall[-1]/wc1])# ; plt.show() ; sys.exit()

axs.set_ylabel(r'$\omega/\Omega_\alpha$',fontsize=20)
axs.set_xlabel(r'$\xi_2$'+'  '+'[%]',fontsize=20)
axs = plt.gca()
axs.ticklabel_format(useOffset=False)
cbar = fig.colorbar(im)
cbar.set_label(r'$\log_{10}(\gamma/\Omega_D)$',fontsize=20)	
#plt.show() ; sys.exit()

## trace through omega to find next curve
intomegas = np.arange(0,35,1)
for i in intomegas:
	if i != 0:
		axs.set_ylim(i-0.5,i+0.5)
	else:
		axs.set_ylim(0,0.5)
	fig.savefig('GAMMA_XI2_{}.png'.format(int(i)),bbox_inches='tight')
	
