
from func_load import *


os.chdir('/storage/space2/phrmsf/lowres_D_He3/')
sims = np.sort([i for i in os.listdir() if 'p_90' in i])
hlabels = [int(i[2:4]) for i in sims]
home = os.getcwd()
# get all sims
for i in range(len(sims)):
	print(sims[i])
	simloc = getSimulation(sims[i])
	ind_lst = list_sdf(simloc)
	
	# quantities and ion species
	name = 'density_energy'
	species = getIonSpecies(sdfread(0))
	quantities = ['Derived_Number_Density_'+i for i in species]
	masses = [getMass(i) for i in species]
	
	try: # load cross cor
		raise SystemError
		crosscorr = read_pkl(name)	
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
		ntot = np.array([masses[i]*narr[i] for i in range(len(species))])
		print(ntot,ntot.shape)
		sys.exit()
		ntot = np.mean([masses[i]*narr[i] for i in range(len(species))],axis=0) # mean across each position		
		plt.imshow(ntot-np.mean(ntot[0:10,:]),**kwargs,cmap='jet')
		plt.xlabel('x',**tnrfont) ; plt.xlabel('t',**tnrfont)
		plt.colorbar()
		plt.savefig(home+'/NumberDensity_{}.png'.format(hlabels[i]))
		plt.clf()
		Ep = read_pkl('Protons_KEdensmatrix')
		dEp = Ep-np.mean(Ep[0:10,:])
		dEp = dEp/np.max(np.abs(dEp))
		# calculate cross correlation
		crosscorr = getCrossCorrelation(ntot,dEp)
	# plot cross correlation
	fig,ax=plotCrossCorrelation(crosscorr,ylabel='t',xlabel='x')
	os.chdir(home)
	fig.savefig('CrossCorrelation_{}.png'.format(hlabels[i]))
	plt.clf()
