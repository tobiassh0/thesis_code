
from func_load import *

# load example sim (densities etc.)
sim_loc = getSimulation('/storage/space2/phrmsf/traceT_highres_0_01')
ind = list_sdf(sim_loc)
d0 = sdfread(0)

# number of compute points
nval = 200000

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
vr = 0.01*v0#*(1/(np.sqrt(2)))
print('vd/v0 :: {},\nvr/v0 :: {}\nu/v_A :: {}'.format(vd/v0,vr/v0,u/vA))

# setup freq range
omegas = wc_maj*np.linspace(0,25,nval)

fig,ax = plt.subplots(figsize=(7,3))
_,k2,_ = coldplasmadispersion(d0,omegas,theta)
posomega, posgamma = growth_rate_man(minions, majions, theta, d0, u, vd, vr, k2, omegas)	
posgamma = np.array(posgamma) #; posgamma[np.isnan(posgamma)] = 0 
knorm = 1/getDebyeLength(d0,'Electrons')
wnorm = wc_maj
ax.plot(posomega/wnorm, np.log10(posgamma/wnorm))
ax.annotate(r'$\theta=$'+str(np.around(theta,1))+r'$^\circ$',xy=(0.1,0.8),xycoords='axes fraction',fontsize=18)
ax.locator_params(axis='y',nbins=4)
ax.locator_params(axis='x',nbins=10)
# 	ax.set_ylabel(r'$\gamma/\Omega_\alpha$',fontsize=18)	
ax.set_ylabel(r'$\gamma/\Omega_\alpha$',fontsize=18)	
ax.set_xlabel(r'$\omega/\Omega_\alpha$',fontsize=18)
#ax.set_ylabel(r'$\gamma/\Omega_\alpha$',fontsize=18)
#fig.add_subplot(111, frameon=False,sharey=ax)
plt.yscale('symlog')
#ax.set_ylim(0,1.4e3)
plt.show()

#os.chdir('/storage/space2/phrmsf/paper/')
#fig.savefig('growthRateMan.png',bbox_inches='tight')
#fig.savefig('growthRateMan.eps',bbox_inches='tight')

