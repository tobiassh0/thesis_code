import numpy as np
from list_new import *
import my_constants as const

fig,ax = plt.subplots(figsize=(8,5))
sims = ['D_89_T_11']#'0_00']#,'D_99_T_01','D_95_T_05','D_89_T_11','D_82_T_18','D_70_T_30','0_50']
xi2 = [0,0.01,0.05,0.11,0.18,0.3,0.5]

for i in range(len(sims)):
	sim_loc = getSimulation('/storage/space2/phrmsf/traceT_'+sims[i])
	d0 = sdfread(0)
	times = read_pkl('times')
	tcD = 2*const.PI/getCyclotronFreq(d0,'Deuterons')
	species = getIonSpecies(d0)
	ENERGYQUANT, labels, fieldmult, fieldquant = getEnergyLabels(sdfread(0),species)
	mean_to = 10
	for i in range(len(ENERGYQUANT)):
		Energyfield = read_pkl(ENERGYQUANT[i])
		meanEnergyfield = np.mean(Energyfield[:mean_to])
		Energyfield = (Energyfield-meanEnergyfield)*fieldmult[i]
		ax.plot(times/tcD,Energyfield)
#	uD = read_pkl('Deuterons_KE')
#	duD = uD - np.mean(uD[:mean_to])
#	uA = read_pkl('Alphas_KE')
#	duA = uA - np.mean(uA[:mean_to])
#	ax.axhline(getMass('Deuterons')/getMass('Alphas'),color='darkgrey',linestyle='--')
#	ax.plot(times/tcD,np.abs(duD/duA))
#	uA = read_pkl('Alphas_KE')
#	duA = uA - np.mean(uA[:mean_to])
#	dt = (times[-1]-times[0])/len(times)
#	sumduA = np.sum(dt*duA)
#	ax.scatter(xi2[i],sumduA,label=xi2[i])

#ax.set_ylabel(r'Total de-energisation'+'  '+'['+r'$Js/m^3$'+']',fontsize=20)
#ax.set_xlabel(r'$\xi_{T}$',fontsize=20)
#fig.savefig('/storage/space2/phrmsf/dump/deenergisation.png',bbox_inches='tight')
#fig.savefig('/storage/space2/phrmsf/dump/deenergisation.eps',bbox_inches='tight')

## normal 
#ax.legend(labels,loc='best')
#ax.set_ylabel(r'$\Delta u$'+'  '+'['+r'$J/m^3$'+']',fontsize=20)
#ax.set_xlabel(r'$t/\tau_{cD}$',fontsize=20)

## EPS 4-page
ax.legend(labels,loc='best',frameon=False)
ax.set_ylabel('Energy  Density'+'  '+'['+r'$J/m^3$'+']',**tnrfont)
ax.set_xlabel(r'Time $t$'+'  '+r'$[\tau{cD}]$',**tnrfont)
fig.savefig('/storage/space2/phrmsf/paper/EPS-4page/eleven_edens.png',bbox_inches='tight',dpi=75)

plt.show()
