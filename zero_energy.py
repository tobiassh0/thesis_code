from lib.list_new import *

fig,ax = plt.subplots(figsize=(8,5))
sims = ['0_00']#,'D_99_T_01','D_95_T_05','D_89_T_11','D_82_T_18','D_70_T_30','0_50']
xi2 = [0,0.01,0.05,0.11,0.18,0.3,0.5]

mean_to = 10
for i in range(len(sims)):
	# plot setup
	fig,ax = plt.subplots(figsize=(10,6),nrows=2,sharex=True)
	fig.subplots_adjust(hspace=0.15)
	
	sim_loc = getSimulation('/storage/space2/phrmsf/traceT_'+sims[i])
	d0 = sdfread(0)
	times = read_pkl('times')
	tcD = 2*const.PI/getCyclotronFreq(d0,'Deuterons')
	species = getIonSpecies(d0)
	ENERGYQUANT, labels, fieldmult, fieldquant = getEnergyLabels(sdfread(0),species)
	print(labels)
	for i in range(len(ENERGYQUANT)):
		Energyfield = read_pkl(ENERGYQUANT[i])/(1000*const.qe)
		if ENERGYQUANT[i] in ['Deuterons_KE','Alphas_KE']:
			tspecies = ENERGYQUANT[i][:-3] # remove '_KE' from name
			pw = getQuantity1d(d0,'Particles_Weight_'+tspecies)
			Nparts = len(pw)*np.mean(pw) # real No. particles = np. simulated * weight per particle
		else:
			Nparts = 1
		meanEnergyfield = np.mean(Energyfield[:mean_to])/Nparts
		Energyfield = ((Energyfield/Nparts-meanEnergyfield)*fieldmult[i])
		if ENERGYQUANT[i] not in ['Deuterons_KE','Alphas_KE']:
			ax[0].plot(times/tcD,Energyfield)
		else:
			print(Energyfield/meanEnergyfield)
			ax[1].plot(times/tcD,Energyfield/meanEnergyfield)
## normal 
ax[0].legend(labels,loc='best')
ax[0].set_ylabel(r'$\Delta u$'+'  '+'['+r'$keV/m^3$'+']',fontsize=20)
ax[0].set_yscale('log')
ax[0].set_ylim(1e12,1e19)
ax[1].legend(labels[3:],loc='best')
ax[1].set_ylabel(r'$\Delta u_\sigma/N_\sigma u_{0\sigma}$',fontsize=20) # +'  '+'['+r'$J/m^3$'+']'
ax[1].set_xlabel(r'$t/\tau_{cD}$',fontsize=20)
ax[1].set_xlim(0,6.1)
fig.savefig('/storage/space2/phrmsf/paper/zero-edens.png')

### EPS 4-page
#ax.legend(labels,loc='best',frameon=False)
#ax.set_ylabel('Energy  Density'+'  '+'['+r'$J/m^3$'+']',**tnrfont)
#ax.set_xlabel(r'Time $t$'+'  '+r'$[\tau{cD}]$',**tnrfont)
#fig.savefig('/storage/space2/phrmsf/paper/EPS-4page/eleven_edens.png',bbox_inches='tight',dpi=75)

plt.show()
