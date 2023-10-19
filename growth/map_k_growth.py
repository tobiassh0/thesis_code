
from func_load import *


#def extractPeaks(data,Nperwcd=1,prominence=0.3):
#	return signal.find_peaks(data,distance=Nperwcd,prominence=prominence)[0] # tune till Nperwcd encapsulates all peaks (visually)

#---------------------------------------------------------------------#
### extract growth rates per omega spacing
#sim_lst = ['traceT_D_89_T_11','traceT_D_99_T_01','traceT_0_00']#'traceT_0_50',
#SimIndex=0
#colors=['g','r','cyan']#'b',
#colors = ['b','r','k']
#shapes = ['o','o','o']
#labels = [r'$11\%$',r'$1\%$',r'$0\%$'] #r'$50\%$',
#fig,axs=plt.subplots(figsize=(10,10/(const.g_ratio)),nrows=2,sharex=True)
#fig.subplots_adjust(hspace=0.1)
#ax = axs.ravel()

##for i in np.arange(0,12,2): ax[0].axvline(i,color='darkgrey',linestyle='--')

#sim_0 = getSimulation('/storage/space2/phrmsf/traceT_0_00')
#wnorm = getCyclotronFreq(sdfread(0),'Deuterons')
#for sim in sim_lst:
#	## define sim
#	sim_loc = getSimulation('/storage/space2/phrmsf/'+sim)
##	wnorm = getEffectiveCyclotronFreq(sdfread(0))
#	omegas, growthRatesMean, growthRatesSTD = map_k_growth(sim_loc,0,25,0.25)#,tstart_frac=0.5,tend_frac=2.0)
##	ax[0].errorbar(omegas/wnorm,growthRatesMean/wnorm,growthRatesSTD/wnorm,fmt='o',color=colors[SimIndex],label=labels[SimIndex])
#	ax[1].scatter(omegas/wnorm,growthRatesMean/wnorm,marker=shapes[SimIndex],facecolors=colors[SimIndex],edgecolors=colors[SimIndex],label=labels[SimIndex])	
#	SimIndex+=1

#ax[1].set_ylabel(r'$\gamma_{s}/\Omega_D$',fontsize=20)	
#ax[1].set_ylim(0,3.5)
##ax[1].set_yscale('symlog')
#ax[1].set_xlim(2,25)
#ax[1].legend(loc='upper left')
#### growth rates
#sim_loc = getSimulation('/storage/space2/phrmsf/traceT_highres_0_01')
#ind = list_sdf(sim_loc)
#nval = 500000

#minions = 'Alphas'
#majions = 'Deuterons'
#elec = 'Electrons'
#theta,_ = getMagneticAngle(sdfread(0))
#wc_maj = getCyclotronFreq(sdfread(0),majions)
#vA = getAlfvenVel(sdfread(0))
#E0 = 3.5E6 * const.eV_to_J
#malpha = const.me*const.me_to_malpha
#v0 = np.sqrt(2*E0/malpha)
#u = np.cos(0.22*const.PI)*v0
#vd = np.sin(0.22*const.PI)*v0
#vr = 0.001*v0#0.001*v0#*(1/(np.sqrt(2)))
#print('vd/v0 :: {},\nvr/v0 :: {}\nu/v_A :: {}'.format(vd/v0,vr/v0,u/vA))
#omegaall = wc_maj*np.linspace(1,25,nval)

#_,kall,_ = coldplasmadispersion(sdfread(0), omegaall)
#posomega, posgamma = growth_rate_man('Alphas', 'Deuterons', 89, sdfread(0), u, vd, vr, kall, omegaall)
#posgamma = np.array(posgamma) #; posgamma[np.isnan(posgamma)] = 0 
#ax[0].plot(posomega/wnorm,posgamma/wnorm,color='k')
#ax[0].set_ylabel(r'$\gamma_{l}/\Omega_D$',fontsize=20)
#ax[0].set_ylim(0,1.4e3)
#ax[1].set_xlabel(r'$\omega/\Omega_D$',fontsize=20)
#plt.show()
#fig.savefig('/storage/space2/phrmsf/paper/Bz_kt_growth_0_25.eps')
#fig.savefig('/storage/space2/phrmsf/paper/Bz_kt_growth_0_25.png')
#---------------------------------------------------------------------#

#---------------------------------------------------------------------#
## multiple growth rates 
sim_lst = ['traceT_D_89_T_11','traceT_D_99_T_01','traceT_0_00']#'traceT_0_50',
colors=['g','r','cyan']#'b',
colors = ['b','r','k']
shapes = ['o','o','o']
labels = [r'$11\%$',r'$1\%$',r'$0\%$'] #r'$50\%$',
fig,axs=plt.subplots(figsize=(8,6),nrows=3,sharex=True)
fig.subplots_adjust(hspace=0.1)
ax = axs.ravel()
times = [[0,0.5],[0.5,2.0]]
wmax = 25

## deuteron harmonics
for a in np.arange(1,axs.shape[0],1):
	for i in np.arange(0,wmax+1,1):
		ax[a].axvline(i,color='darkgrey',linestyle='--')

## theory
minions = 'Alphas'
majions = 'Deuterons'
elec = 'Electrons'
simloc = getSimulation('/storage/space2/phrmsf/traceT_0_00')
d0 = sdfread(0)
theta,_ = getMagneticAngle(d0)
wc_maj = getCyclotronFreq(d0,majions)
wc_min = getCyclotronFreq(d0,minions)
wnorm = wc_min
vA = getAlfvenVel(d0)
E0 = 3.5E6 * const.qe
malpha = const.me*const.me_to_malpha
v0 = np.sqrt(2*E0/malpha)
vperp0 = np.cos(0.22*const.PI)*v0
vpara0 = np.sin(0.22*const.PI)*v0
vr = 0.001*v0#0.01*v0#*(1/(np.sqrt(2)))
print('vpara0/v0 :: {},\nvr/v0 :: {}\nvperp0/v_A :: {}'.format(vpara0/v0,vr/v0,vperp0/vA))
nval = 200000
omegaall = wc_maj*np.linspace(1,25,nval)

_,kall,_ = coldplasmadispersion(d0, omegaall)
posomega, posgamma = growth_rate_man('Alphas', 'Deuterons', 89, d0, vperp0, vpara0, vr, kall, omegaall)
posgamma = np.array(posgamma) #; posgamma[np.isnan(posgamma)] = 0 
ax[0].plot(posomega/wnorm,posgamma/wnorm,color='k')

## empirical growths
t=1
for ttimes in times:
	SimIndex=0
	for sim in sim_lst:
		## define sim
		sim_loc = getSimulation('/storage/space2/phrmsf/'+sim)
		omegas, growthRatesMean, growthRatesSTD = map_k_growth(sim_loc,'Deuterons',0,wmax,0.25,tstart_frac=ttimes[0],tend_frac=ttimes[1])
#		ax[t].scatter(omegas/wnorm,growthRatesMean/wnorm,marker='o',facecolors=colors[SimIndex])
		thresh = growthRatesMean > 0
		growthRatesMean = growthRatesMean[thresh]
		omegas = omegas[thresh]
		ax[t].plot(omegas/wnorm,growthRatesMean/wnorm,'-o',color=colors[SimIndex])
		SimIndex+=1
	ax[t].set_ylim(0,3.5)
	ax[t].set_xlim(0,wmax)
	ax[t].set_ylabel(r'$\gamma_s/\Omega_\alpha$',**tnrfont)
	t+=1

# time range annotations
ax[0].set_ylabel(r'$\gamma_l/\Omega_\alpha$',**tnrfont)
ax[1].text(1.,2.8,r'$0<t/\tau_{cD}<0.5$',color='black',fontsize=16,bbox=dict(facecolor='white',edgecolor='black',boxstyle='square,pad=0.25'))
ax[2].text(1.,2.8,r'$0.5<t/\tau_{cD}<2.0$',color='black',fontsize=16,bbox=dict(facecolor='white',edgecolor='black',boxstyle='square,pad=0.25'))
# x-label
ax[2].set_xlabel(r'$\omega/\Omega_\alpha$',**tnrfont)
#fig.savefig('/storage/space2/phrmsf/paper/Bz_kt_growth_twoTimes.png')
#fig.savefig('/home/space/phrmsf/Documents/thesis_code/Bz_kt_growth_twoTimes.png')
plt.show()	

#---------------------------------------------------------------------#

#---------------------------------------------------------------------#
##def PITENSOR(file0,v0,kall,omegaall,theta=90):
##	theta = theta*(const.PI/180) # radians
##	PIxx = np.zeros(len(omegaall),complex) ; PIxy = np.zeros(len(omegaall),complex) ; PIyy = np.zeros(len(omegaall),complex)
##	wci = getEffectiveCyclotronFreq(file0)

##	for i in range(0, omegaall.shape[0]):
##		l = round(omegaall[i]/wci) #l closest to the omega
##		k = kall[i]
##		kpara = kall[i]*np.cos(theta)
##		kperp = kall[i]*np.sin(theta)
##		
##		za = kperp*v0/wci
##		zarr = np.linspace(0,2*za,10000)
##		dzarr = (zarr[-1]-zarr[0])/len(zarr)

##		# try and except clause to remove l=0 division
##		try:
##			PIxx[i] = spec.jv(2*l,2*za)
##		except:
##			None
##		try:
##			PIxy[i] = (-1j*za/l)*(spec.jvp(2*l,2*za))
##		except:
##			None
##		try:
##			PIyy[i] = (1-(za/l)**2)*spec.jv(2*l,2*za)+(za/(2*l**2))*np.sum(dzarr*spec.jv(2*l,zarr))
##		except:
##			None

##	return PIxx, PIxy, PIyy

##def CHI0CALC(PIxx,PIxy,PIyy,n0,kall,omegaall,M,L,Z1,Z2,Z3,species=['Deuterons','Tritons','Alphas'],B=2.1):
##	meff = (getMass(species[0])/Z1)*(1-M*Z2-L*Z3) + getMass(species[1])*M + getMass(species[2])*L
##	print(PIxx.shape,PIxy.shape,PIyy.shape,meff.shape,omegaall.shape,kall.shape)
##	Chi02 = PIxx + (const.qe*const.mu0*n0/(kall))*(((-2*1j)*omegaall/B)*PIxy + (1/meff)*PIyy)
##	return np.sqrt(Chi02)


##sim_loc = getSimulation('/storage/space2/phrmsf/traceT_D_99_T_01')
##nval = 100
##Mu = np.linspace(0,1,100)
##Lambda = np.linspace(1e-5,1e-2,100)
##Z1 = 1 ; Z2 = 1 ; Z3 = 2
##CHI0 = np.zeros((len(Mu),len(Lambda),10000))
##species = ['Deuterons','Tritons','Alphas']
##B0 = 2.1
##E0 = 3.5E6 * const.qe
##v0 = np.sqrt(2*E0/getMass('Alphas'))
##n0 = getMeanquantity(sdfread(0),'Derived_Number_Density_Electrons')

##for m in range(len(Mu)):
##	for l in range(len(Lambda)):
##		meff = (getMass(species[0])/Z1)*(1-Mu[m]*Z2-Lambda[l]*Z3) + getMass(species[1])*Mu[m] + getMass(species[2])*Lambda[l]
##		wci = const.qe*B0/meff
##		omegas = wci*np.linspace(0,25,10*nval)
##		_,k2,_ = coldplasmadispersion(sdfread(0),omegas)
##		PIxx, PIxy, PIyy = PITENSOR(sdfread(0),v0,k2,omegas,theta=89)
##		CHI0[m,l,:] = CHI0CALC(PIxx,PIxy,PIyy,n0,k2,omegas,Mu[m],Lambda[l],Z1=1,Z2=1,Z3=2,species=['Deuterons','Tritons','Alphas'],B=2.1)
##dumpfiles(CHI0,'CHI0_Mu_Lambda')
##plt.imshow(CHI0,**kwargs) ; plt.show()




