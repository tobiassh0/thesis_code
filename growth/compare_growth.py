
from func_load import *


kwargs={'interpolation':'nearest','origin':'lower','aspect':'auto'}


def coldplasma_eff(file0,omegas,theta=None):
	if not theta: 
		theta, _ = getMagneticAngle(file0) # assuming B directed in z-x plane
	else:
		theta = theta*const.PI/180
	sin = np.sin(theta) ; cos = np.cos(theta)
	print(theta, sin, cos)
	
	wci = getEffectiveCyclotronFreq(file0)
	wpi = getEffectivePlasmaFreq(file0)
	wpe = getPlasmaFreq(file0,'Electrons')
	wce = getCyclotronFreq(file0,'Electrons',getChargeNum('Electrons'))
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
def growth_rate_Eff(minions, theta, b, u, vd, vr, kall, omegaall):

	wcycb = getCyclotronFreq(b,minions) # cyc freq for beam ions
	wcyci = getEffectiveCyclotronFreq(b) # cyc freq for bulk ions

	wpb = getPlasmaFreq(b, minions)
	wpi = getEffectivePlasmaFreq(b)
	va = getAlfvenVel(b)

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
		if (gammas[i] >= 0):
			posomega.append(omegaall[i])
			posgamma.append(gammas[i])

	return np.array(posomega,dtype='float'), np.array(posgamma,dtype='float')


theta = 89
B0 = 2.1
sim_loc = getSimulation('/storage/space2/phrmsf/traceT_D_89_T_11')
file0 = sdfread(0)
times = read_pkl('times')
wcycD = const.qe*B0/getMass('Deuterons')
nval = 200000
omegaall = wcycD * np.linspace(0,50,nval)

E0 = 3.5E6 * const.eV_to_J
malpha = const.me*const.me_to_malpha
v0 = np.sqrt(2*E0/malpha)
u = np.cos(0.22*const.PI)*v0
vd = np.sin(0.22*const.PI)*v0
vr = 0.001*v0

colors = ['r','b','k']
knorm = 1/getDebyeLength(file0,'Electrons') # normalises to units of radians (# of 2pi)
dx = getdxyz(file0) ; dt = times[-1]/len(times)
klim = 0.5*2*const.PI/dx ; wlim = 0.5*2*const.PI/dt
FT2d = read_pkl('FT_2d_Magnetic_Field_Bz')

### compare growth rates
fig,axs=plt.subplots(ncols=4,figsize=(12,6),sharey=True)
fig.subplots_adjust(wspace=0.05)
ax=axs.ravel()
ax[0].imshow(np.log10(FT2d),**kwargs,cmap='magma',extent=[0,klim/knorm,0,wlim/wcycD])
_,k2eff,_=coldplasma_eff(file0, omegaall, theta)
posomegaEff,posgammaEff = growth_rate_Eff('Alphas', theta, file0, u, vd, vr, k2eff, omegaall)
_,k2single,_=coldplasmadispersion(file0, 'Deuterons', '', omegaall, theta)
posomegaSingle, posgammaSingle = growth_rate_man('Alphas', 'Deuterons', theta, file0, u, vd, vr, k2single, omegaall)
_,k2two,_=coldplasmadispersion(file0, 'Deuterons', 'Tritons', omegaall, theta)
posomegaTwo, posgammaTwo = growth_rate_man('Alphas', 'Deuterons', theta, file0, u, vd, vr, k2two, omegaall)

ax[0].plot(k2eff/knorm,omegaall/wcycD,color=colors[0],label='eff')
ax[0].plot(k2single/knorm,omegaall/wcycD,color=colors[1],label='one')
ax[0].plot(k2two/knorm,omegaall/wcycD,color=colors[2],label='two')

ax[1].plot(posgammaSingle/wcycD,posomegaSingle/wcycD,color=colors[1])
ax[2].plot(posgammaEff/wcycD,posomegaEff/wcycD,color=colors[0])
ax[3].plot(posgammaTwo/wcycD,posomegaTwo/wcycD,color=colors[2])

#ax[0].legend(loc='best')
ax[0].set_xlim(0,0.06) ; ax[0].set_ylim(0,32)
ax[1].set_ylim(0,32)
#ax[1].set_xscale('symlog')
ax[0].set_xlabel(r'$k\lambda_{De}$',fontsize=20)
fig.text(0.075, 0.5, r'$\omega/\Omega_D$', va='center', rotation='vertical',fontsize=20)
#ax[0].set_ylabel(r'$\omega/\Omega_D$',fontsize=20)
text = ['','One ion','Effective','Two ions']

ax[0].set_xticks([0,0.015,0.03,0.045])
ax[0].set_xticklabels(['0','0.015','0.03','0.045'])
for i in [1,2,3]:
	ax[i].set_xlim(0,2500)
	ax[i].set_xlabel(r'$\gamma/\Omega_D$',fontsize=20)
	ax[i].annotate(text[i],(0.05,0.95),xycoords='axes fraction')
	ax[i].set_xticks([0,750,1500,2250])
	ax[i].set_xticklabels(['0','750','1500','2250'])
	ax[i].set_yticklabels([])
ax[0].set_yticklabels([0,5,10,15,20,25,30])	

fig.savefig('/storage/space2/phrmsf/dump/CompareGrowth.png',bbox_inches='tight')
plt.clf()

############################################################################

#### compare quality of each prediction 
fig,axs=plt.subplots(nrows=2,sharex=True)
fig.subplots_adjust(hspace=0.075)
ax=axs.ravel()
ax[0].imshow(np.log10(FT2d),**kwargs,cmap='magma',extent=[0,klim/knorm,0,wlim/wcycD])
## extract strongest FAW k per w
#karr = np.linspace(0,klim,FT2d.shape[1])
#warr = (FT2d.argmax(axis=0))*wlim/FT2d.shape[0] 
warr = np.linspace(0,wlim,FT2d.shape[0])
karr = (FT2d.argmax(axis=1))*klim/FT2d.shape[1]
print(warr.shape,karr.shape)
#ax[0].scatter(karr/knorm,warr/wcycD,color='k',marker='x')

## find difference between predicted k (from each model) and strongest k
# limit comparison up to k/knorm < 0.15
thresh = (warr/wcycD<32)
warr = warr[thresh]
karr = karr[thresh]
omegas = np.linspace(0,warr[-1],len(warr))
_,k2eff,_=coldplasma_eff(file0, omegas, theta)
_,k2single,_=coldplasmadispersion(file0, 'Deuterons', '', omegas, theta)
_,k2two,_=coldplasmadispersion(file0, 'Deuterons', 'Tritons', omegas, theta)
ax[0].plot(k2eff/knorm,omegas/wcycD,color=colors[0],label='eff')
ax[0].plot(k2single/knorm,omegas/wcycD,color=colors[1],label='one')
ax[0].plot(k2two/knorm,omegas/wcycD,color=colors[2],label='two')

erreff = np.zeros(len(omegas))
errsingle = np.zeros(len(omegas))
errtwo = np.zeros(len(omegas))
for i in range(len(omegas)):
	erreff[i] = np.sqrt((((karr[i]-k2eff[i])/knorm)**2) + (((warr[i]-omegas[i])/wcycD)**2))
	errsingle[i] = np.sqrt((((karr[i]-k2single[i])/knorm)**2) + (((warr[i]-omegas[i])/wcycD)**2))
	errtwo[i] = np.sqrt((((karr[i]-k2two[i])/knorm)**2) + (((warr[i]-omegas[i])/wcycD)**2))
ax[1].scatter(k2eff/knorm,erreff,color=colors[0],marker='x',label='Eff')
ax[1].scatter(k2single/knorm,errsingle,color=colors[1],marker='x',label='One')
ax[1].scatter(k2two/knorm,errtwo,color=colors[2],marker='x',label='Two')

meanerreff = np.nanmean(erreff) ; stderreff = np.nanstd(erreff)
meanerrsingle = np.nanmean(errsingle) ; stderrsingle = np.nanstd(errsingle)
meanerrtwo = np.nanmean(errtwo) ; stderrtwo = np.nanstd(errtwo)
print(str(meanerreff)+' +- '+str(stderreff))
print(str(meanerrsingle)+' +- '+str(stderrsingle))
print(str(meanerrtwo)+' +- '+str(stderrtwo))

ax[1].legend(loc='best')
ax[0].set_ylabel(r'$\omega/\Omega_D$',fontsize=20)
ax[0].set_ylim(0,32)
ax[1].set_ylabel(r'$d_{k\omega}$',fontsize=20)
ax[1].set_xlabel(r'$k\lambda_{De}$',fontsize=20)
ax[1].set_xlim(0,karr[-1]/knorm)
ax[1].set_yscale('log')

fig.savefig('/storage/space2/phrmsf/dump/CompareColdPlasmaDispersion.png')
plt.show()



