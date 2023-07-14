from list_new import *

def Endict(field):
	Energydict = {
	'Magnetic_Field_Bz': 'Bzenergy',
	'Magnetic_Field_By': 'Byenergy',
	'Magnetic_Field_Bx': 'Bxenergy',
	'Electric_Field_Ez': 'Ezenergy',
	'Electric_Field_Ey': 'Eyenergy',
	'Electric_Field_Ex': 'Exenergy',
	}
	EnergyLabeldict = {
	'Magnetic_Field_Bz': r'$\langle\Delta B_z^2\rangle/2\mu_0$',
	'Magnetic_Field_By': r'$\langle B_y^2\rangle/2\mu_0$',
	'Magnetic_Field_Bx': r'$\langle B_x^2\rangle/2\mu_0$',
	'Electric_Field_Ez': r'$\langle E_z^2\rangle\epsilon_0/2$',
	'Electric_Field_Ey': r'$\langle E_y^2\rangle\epsilon_0/2$',
	'Electric_Field_Ex': r'$\langle E_x^2\rangle\epsilon_0/2$',
	}
	Energydict = {
	'Magnetic_Field_Bz': const.mult_magn_field
	'Magnetic_Field_By': const.mult_magn_field,
	'Magnetic_Field_Bx': const.mult_magn_field,
	'Electric_Field_Ez': const.mult_elec_field,
	'Electric_Field_Ey': const.mult_elec_field,
	'Electric_Field_Ex': const.mult_elec_field,
	}

	return Energydict.get(field), EnergyLabeldict.get(field), Energymult.get(field)
	
def getFieldEnergy(IFN):
	index, F0, name = IFN
	return np.mean((getQuantity1d(sdfread(index),name)-F0)**2)

def getParticleEnergy(IFN):
	index, F0, name = IFN
	ek = getQuantity1d(sdfread(index),'Derived_Average_Particle_Energy_'+name)
	dens = getQuantity1d(sdfread(index),'Derived_Number_Density_'+name)
	return np.mean(ek*dens)
	
def getEnergy(ARR):
	index, F0, names = ARR
	Energy = []
	for name in names:
		if 'Field' in name:
			return np.mean((getQuantity1d(sdfread(index),name)-F0)**2)
		else:
			ek = getQuantity1d(sdfread(index),'Derived_Average_Particle_Energy_'+name)
			dens = getQuantity1d(sdfread(index),'Derived_Number_Density_'+name)
			return np.mean(ek*dens)

# Home function which can be called, will call on getEnergies(). Generates KE_Dens and field E_Dens for species 1, 2, 3 and e-
def ENERGIES(sim_loc,frac=1,mean_to=10):
	width = 10
	height = 6
	#height = width*(1/const.g_ratio)
	fig,ax=plt.subplots(figsize=(width,height))
	
	## check if field values in normal dump 
	try:
		d = sdfread(1) # hardcoded for now, will default if cant load
		keys = getKeys(d)
		fieldquant = []
		for key in keys:
			if 'Field' in key:
				fieldquant.append(key) # append if is a field value
	except:
		fieldquant=['Magnetic_Field_Bz','Magnetic_Field_By','Electrical_Field_Ex'] # default fields

	## append to energy_quant array
	energy_quant=[] ; names=[] ; Energy_mult=[]
	for field in fieldquant:
		fieldlabel, labelval, Emult = Endict(field)
		energy_quant.append(fieldlabel)
		names.append(labelval)
		Energy_mult.append(Emult)
	
	species = getIonSpecies(sdfread(1))
	maj_species, maj2_species, min_species = species
	species = np.append(species,'Electrons')
	ionquant = [] ; ionlabels = []

	for spec in species:
		if spec != '': 
			ionquant.append(spec+'_KE')
			ionlabels.append(getIonlabel(spec))
	
	for k in range(len(ionquant)):
		fieldquant.append(species[k])
		energy_quant.append(ionquant[k])
		names.append(ionlabels[k])

	## do when no minority ring-beam
	if maj2_species == '': 
		if min_species == '': # single ion
			Single = True
			Double = Triple = False
		else: # 2 ions:
			Double = True
			Single = Triple = False
	else: # 3 ions
		Triple = True
		Single = Double = False

	Zmaj1 = getChargeNum(maj_species)
	Zmaj2 = getChargeNum(maj2_species)
	Zmin = getChargeNum(min_species)

	index_list = list_sdf(sim_loc)
#	os.chdir(sim_loc)
	times = read_pkl('times')
	n = (len(index_list)//frac)

	F0 = np.zeros(len(energy_quant))
	for i in range(len(energy_quant)):
		if energy_quant[i] == 'Magnetic_Field_Bz':
			F0[i] = getQuantity1d(sdfread(0), 'Magnetic_Field_Bz')

	Energies = np.zeros((len(fieldquant),n))	
	for t in range(0,n):
		if t%(n//20)==0: print(str(round(100*t/n))+' %') # print every 5%
		d = sdfread(t)
		for s in range(len(fieldquant))
			if 'Field' in fieldquant[s]:
				Energies[s,t] = np.mean((getQuantity1d(d,fieldquant[s])-F0[s])**2)
			else:
				Energies[s,t] = getTotalKineticEnergyDen(d,fieldquant[s])

	## energy_quant and fieldquant should be the same length 
	for s in range(len(energy_quant)):
		dumpfiles(Energies[s,:],energy_quant[s])


############################################
#	for i in range(len(fieldquant)): # loops through fields and particles
#		Energy = np.zeros(n)
#		if 'Bz' in fieldquant[i]: F0 = getQuantity1d(sdfread(0),fieldquant[i])
#		else: F0 = 0
#		tindex_list = [[j, F0, fieldquant[i]] for j in index_list]
##		print(tindex_list)
#		if 'Field' in fieldquant[i]:
#			func = getFieldEnergy
#		else:
#			func = getParticleEnergy
#		pool = mp.Pool(mp.cpu_count())
#		Energy = np.array(pool.map_async(func,tindex_list).get(99999))
#		pool.close()
#		dumpfiles(Energy,energy_quant[i]+'_para')
##		Energy = Energy - np.mean(Energy[:mean_to])
##		plt.plot(Energy) ; plt.show()	
##########################################################


#########################################################
#########################################################
import time

sim_loc = getSimulation('/storage/space2/phrmsf/traceT_0_01')
tic = time.time()
eint = energies(sim_loc,min_species='Alphas',maj_species='Deuterons',maj2_species='Tritons',frac=1,plot=False,integ=False)
toc = time.time()
tloop = (toc-tic)

tic = time.time()
ENERGIES(sim_loc)
toc = time.time()
tpara = (toc-tic)

dumpfiles([tpara,tloop],'ENERGYTEST_TIMES')
print('PARALLEL :: {}\nLOOP :: {}'.format(tpara,tloop))



