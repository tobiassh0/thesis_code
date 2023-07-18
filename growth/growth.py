
from func_load import *


sim_loc = getSimulation('/storage/space2/phrmsf/traceT_highres_0_01')
ind = list_sdf(sim_loc)
nval = 200000

minions = 'Alphas'
majions = 'Deuterons'
elec = 'Electrons'
theta,_ = getMagneticAngle(sdfread(0))
wc_maj = getCyclotronFreq(sdfread(0),majions)
vA = getAlfvenVel(sdfread(0))
E0 = 3.5E6 * const.eV_to_J
malpha = const.me*const.me_to_malpha
v0 = np.sqrt(2*E0/malpha)
u = np.cos(0.22*const.PI)*v0
vd = np.sin(0.22*const.PI)*v0
vr = 0.06*v0#0.001*v0#*(1/(np.sqrt(2)))
print('vd/v0 :: {},\nvr/v0 :: {}\nu/v_A :: {}'.format(vd/v0,vr/v0,u/vA))
omegas = wc_maj*np.linspace(0,25,nval)
thetas = [89]


#fig, axs = plt.subplots(len(thetas), sharex=True)
#plt.subplots_adjust(wspace=1.)
#for i in range(len(thetas)):
#	_,k2,_ = coldplasmadispersion(sdfread(0), majions, elec, getChargeNum(majions), getChargeNum(elec), omegas, thetas[i])
#	posomega, posgamma = growth_rate_man(minions, majions, thetas[i], sdfread(0), v0, vd, vr, k2, omegas)
#	
#	ax = axs[i]
#	ax.plot(posomega, posgamma)
##	axs[i].set_ylim(0,1.2)
#	ax.annotate(r'$\theta=$'+str(np.around(thetas[i],1))+r'$^\circ$',xy=(0.015,0.75),xycoords='axes fraction')
#	ax.locator_params(axis='y',nbins=4)
#	ax.locator_params(axis='x',nbins=10)
#	# if thetas[i] == thetas[len(thetas)//2]: # middle value
#	# 	ax.set_ylabel(r'$\gamma/\Omega_\alpha$',fontsize=18)	

#ax.set_ylabel(r'$\gamma/\Omega_\alpha$',fontsize=18)	
##fig.add_subplot(111, frameon=False,sharey=ax)
#plt.xlabel(r'$\omega/\Omega_\alpha$',fontsize=18)
##plt.ylabel(r'$\gamma/\Omega_\alpha$',fontsize=18)
#plt.show()
##fig.savefig('Growth_rates_JET26148.jpg',bbox_inches='tight')

fig,ax = plt.subplots(figsize=(7,3))
_,k2,_ = coldplasmadispersion(sdfread(0),majions,'',omegas)
posomega, posgamma = growth_rate_man(minions, majions, thetas[0], sdfread(0), u, vd, vr, k2, omegas)	
posgamma = np.array(posgamma) #; posgamma[np.isnan(posgamma)] = 0 
knorm = 1/getDebyeLength(sdfread(0),'Electrons')
wnorm = wc_maj
ax.plot(posomega/wnorm, posgamma/wnorm)
ax.annotate(r'$\theta=$'+str(np.around(thetas[0],1))+r'$^\circ$',xy=(0.1,0.8),xycoords='axes fraction',fontsize=18)
ax.locator_params(axis='y',nbins=4)
ax.locator_params(axis='x',nbins=10)
# if thetas[i] == thetas[len(thetas)//2]: # middle value
# 	ax.set_ylabel(r'$\gamma/\Omega_\alpha$',fontsize=18)	
ax.set_ylabel(r'$\gamma/\Omega_\alpha$',fontsize=18)	
ax.set_xlabel(r'$\omega/\Omega_\alpha$',fontsize=18)
#ax.set_ylabel(r'$\gamma/\Omega_\alpha$',fontsize=18)
#fig.add_subplot(111, frameon=False,sharey=ax)
#plt.yscale('symlog')
#ax.set_ylim(0,1.4e3)
plt.show()

#os.chdir('/storage/space2/phrmsf/paper/')
#fig.savefig('growthRateMan.png',bbox_inches='tight')
#fig.savefig('growthRateMan.eps',bbox_inches='tight')





