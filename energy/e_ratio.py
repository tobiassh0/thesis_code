#import matplotlib.pyplot as plt 
from func_load import *


## loop through simulations and species (num of panels)
figpart, axpart = plt.subplots(nrows=3,figsize=(8,6),sharex=True) # assume 3 ion species # (6,10)
figpart.subplots_adjust(hspace=0.17)
## fields
figfield, axfield = plt.subplots(nrows=1,ncols=3,figsize=(10,6)) # assume 3 ion species # nrows=2,ncols=3figsize=(18,8)
figfield.subplots_adjust(wspace=0.,hspace=0.1)

mean_to = 10
#sim_lst = ['traceT_0_50','traceT_D_70_T_30','traceT_D_82_T_18','traceT_D_89_T_11','traceT_D_95_T_05','traceT_D_99_T_01','traceT_0_00']
sim_lst = ['traceT_0_50','traceT_D_89_T_11','traceT_D_99_T_01','traceT_0_00']
fields = ['Electric_Field_Ex','Magnetic_Field_By','Magnetic_Field_Bz']
fieldnames = [] ; fieldmult = [] ; labels = []
for field in fields:
	a,b,c = Endict(field)
	fieldnames.append(a)
	labels.append(b)
	fieldmult.append(c)

#colors = ['blue','pink','orange','g','k','r','cyan']
colors = ['blue','g','r','cyan']

c = 0
N = 300
#alphas = [0.1,0.33,0.66,0.99]
#colors = ['lightblue','lightskyblue','dodgerblue','royalblue']
growth=[]
for sim in sim_lst:
	print(sim)
	sim_loc = getSimulation('/storage/space2/phrmsf/'+sim)
	times = read_pkl('times')	
	d0 = sdfread(0)
	nx = len(getQuantity1d(d0,'Derived_Number_Density'))
	n0 = getQuantity1d(d0,'Derived_Number_Density_Electrons')
	wcD = getCyclotronFreq(d0,'Deuterons')
	tcD = 2*const.PI/wcD
	dt = (times[-1]-times[0])/len(times)
	dt_p = dt/tcD
	print('################## :: ',dt,tcD,dt_p)
	species = getIonSpecies(d0)
	Nparts = np.zeros(len(species))
	## field & species energies
	for i in range(0,3):
		# fields
		Energyfield = read_pkl(fieldnames[i])
		if fieldnames[i] == 'Bzenergy':
			meanEnergyfield = np.mean(Energyfield[:mean_to])
			fitting = True
		elif fieldnames[i] == 'Alphas':
			## running mean
			Energyfield = np.convolve(Energyfield,np.ones(N)/N,mode='valid')
			timesplot = np.linspace(0,max(times),len(Energyfield))
			axfield[i].plot(timesplot/tcD,Energyfield,color=colors[c])
#			## gradient of field ===> d(Delta u)/dt
#			gradDu = np.convolve(np.gradient((Energyfield),dt),np.ones(N)/N,mode='valid') # np.convolve
#			timesplot = np.linspace(0,max(times),len(gradDu))
#			axfield[1,i].plot(timesplot/tcD,np.log(gradDu),color=colors[c])
		else:
			meanEnergyfield = 0
		Energyfield = (Energyfield-meanEnergyfield)*fieldmult[i]
		## running mean
		Energyfield = np.convolve(Energyfield,np.ones(N)/N,mode='valid')
		timesplot = np.linspace(0,max(times),len(Energyfield))
		axfield[i].plot(timesplot/tcD,np.log(Energyfield),color=colors[c])
		## gradient of field ===> d(Delta u)/dt
#		gradDu = np.convolve(np.gradient((Energyfield),dt),np.ones(N)/N,mode='valid') # np.convolve
#		timesplot = np.linspace(0,max(times),len(gradDu))
#		axfield[1,i].plot(timesplot/tcD,np.log(gradDu),color=colors[c])
#		axfield[1,i].set_xlim(0,2.5)
#		axfield[1,i].set_ylim(13,23)
		axfield[i].annotate(labels[i],xy=[0.1,0.15],xycoords='axes fraction',fontsize=20)#,bbox=dict(fc="white"))
#		axfield[1,i].annotate(labels[i],xy=[0.025,0.85],xycoords='axes fraction',fontsize=20)#,bbox=dict(fc="white"))
		axfield[i].set_ylim(-6,6) # axfield[i].set_ylim(-0.1,265)
		axfield[i].set_xticklabels([0,1,2,3,4,5,6])
#		axfield[1,i].set_xticklabels([0,0.5,1.0,1.5,2.0])
		if i == 0:
			axfield[i].set_ylabel(r'$\ln[\Delta u]$'+' '+'['+r'$Jm^{-3}$'+']',fontsize=20)
#			axfield[i].set_ylabel(r'Energy Density'+'  '+'['+r'$Jm^{-3}$'+']',**tnrfont)
		else: # i > 0
			axfield[i].set_yticklabels([])
#			axfield[1,i].set_yticklabels([])
		## print growth rates for fields and d(Du)/dt
		thresh = (timesplot > 0) & (timesplot < 0.5*tcD)
		timefit = timesplot[thresh]
		datafit = Energyfield[thresh]
		popt, pcov = curve_fit(lambda t,g,off: np.exp(g*t)+off,timefit,datafit,maxfev=5000) # exponential fitting
		g,off = popt
		growth.append([sim,fieldnames[i],g,off])

		# particles
		try:
			pw = getQuantity1d(d0,'Particles_Weight_'+species[i])
			Nparts[i] = 1#len(pw)*np.mean(pw) # real No. particles = np. simulated * weight per particle
			nspec = getMeanquantity(d0,'Derived_Number_Density_'+species[i])
		except:
			nspec = 1
			continue # should be 0 for species that arent present
		Energypart = read_pkl(species[i]+'_KE')/(nspec*1000*const.qe) # energy density 
		meanEnergypart = np.mean(Energypart[:mean_to])
		print(Energypart,meanEnergypart)
		Energypart = np.convolve(Energypart,np.ones(N)/N,mode='valid')
		timespart = np.linspace(0,max(times),len(Energypart))
		thresh = timespart/tcD < 6.1
		print(species[i],Nparts[i])
		axpart[i].plot(timespart[thresh]/tcD,(Energypart[thresh]-meanEnergypart)/Nparts[i],color=colors[c])#alpha=alphas[c]) #
		axpart[i].set_xlim(0,6.1)
		if species[i] in ['Deuterons', 'Tritons']: 
			axpart[i].set_ylim(-0.02,0.3)
		else:
			axpart[i].set_ylim(-300,10)
	c+=1


#growth=np.array(growth).reshape(-1,4)
#print(growth,growth.shape)
#print('wcD :: ',wcD)
#totgrowth=0
#for i in range(growth.shape[0]):
#	if growth[i,1]=='Bzenergy':
#		totgrowth+=float(growth[i,2])
#	else:
#		None
#meangrowth = totgrowth/int(growth.shape[0])
#print('meangrowth [wcD] :: ',2*const.PI*meangrowth/wcD)

labels = [r'$50\%$',r'$30\%$',r'$18\%$',r'$11\%$',r'$5\%$',r'$1\%$',r'$0\%$']
labels = [r'$50\%$',r'$11\%$',r'$1\%$',r'$0\%$']
axfield[1].legend(labels,loc='best',frameon=False)
axpart[2].legend(labels,loc='best')#,frameon=False)
labels = [r'Deuterons',r'Tritons',r'Alphas']
labels = [r'$\Delta E_D$',r'$\Delta E_{T}$',r'$\Delta E_\alpha$']

for i in range(len(axfield)):
	axfield[i].set_xlabel(r'$t/\tau_{cD}$',fontsize=20)
#	axfield[i].set_xlabel(r'Time $t$'+'  '+r'$[\tau_{cD}]$',**tnrfont)
for j in range(len(axpart)): 
	axpart[j].set_ylabel(labels[j]+' '+'['+r'$keV$'+']',**tnrfont)
#	axpart[j].annotate(labels[j],xy=(0.75,0.8),xytext=None,xycoords='axes fraction',ha='left',**tnrfont)


## 4-page EPS
#axpart[1].set_ylabel('Energy Density per-particle'+'  '+r'$[keV$'+ ' '+r'$m^{-3}]$',**tnrfont)
#axpart[-1].set_xlabel(r'Time $t$'+'  '+r'$[\tau_{cD}]$',**tnrfont)
#figfield.savefig('/storage/space2/phrmsf/paper/EPS-4page/fieldEnergies.png',bbox_inches='tight',dpi=75)
#figpart.savefig('/storage/space2/phrmsf/paper/EPS-4page/Energy_Dens_pp_v2.png',bbox_inches='tight',dpi=75)

## normal
axpart[-1].set_xlabel(r'$t/\tau_{cD}$',fontsize=20)
#figfield.savefig('/storage/space2/phrmsf/paper/fieldEnergies.eps',bbox_inches='tight')
figpart.savefig('/storage/space2/phrmsf/paper/Energy_Dens_pp_v2.eps')
figpart.savefig('/storage/space2/phrmsf/paper/Energy_Dens_pp_v2.png')

#plt.show()











