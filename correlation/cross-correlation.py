
from func_load import *

# get sim
simloc = getSimulation('/storage/space2/phrmsf/lowres_D_He3/0_34_p_90')
ind_lst = list_sdf(simloc)

# quantities and ion species
name = 'density_energy'
species = getIonSpecies(sdfread(0))
quantities = ['Derived_Number_Density_'+i for i in species]
masses = [getMass(i) for i in species]

try: # load cross cor
	crosscor = read_pkl(name)	
except: # calculate cross cor
	# get number density and energy
	narr = []
	try:
		for j in range(len(species)):
			narr.append(load_batch_fieldmatrix([],quantities[j],para=False))
	except:
		get_batch_fieldmatrix(ind_lst,quantities=quantities,load=False)
		for j in range(len(species)):
			narr.append(load_batch_fieldmatrix([],quantities[j],para=False))
	
	# weighted number density according to mass
	ntot = np.mean([masses[i]*narr[i] for i in range(len(species))],axis=0) # mean across each position
	
Ep = read_pkl('Protons_KEdensmatrix')
dEp = Ep-np.mean(Ep[0:10,:])

# calculate cross correlation
crosscorr = getCrossCorrelation(ntot,dEp)
plotCrossCorrelation(crosscor,ylabel='t',xlabel='x')