
from func_load import *

## get the cross correlation between two matrices
def getCrossCorrelationMat(mat1,mat2,name='density_energy',axis=0):
	if mat1.shape != mat2.shape:
		print('# ERROR # :: Cant calculate cross-correlation of unequal arrays')
		sys.exit()
	else:
		# setup cross correlation matrix
		crosscor = np.zeros(mat1.shape)
		# might not match up
		l = mat1.shape[axis]
		print(l)
		if axis == 0:
			for i in range(l):
				print(i/l*100)
				# normalise both matrices according to https://en.wikipedia.org/wiki/Cross-correlation#Zero-normalized_cross-correlation_(ZNCC)
				normat1 = (mat1[i,:]-np.mean(mat1[i,:]))/(np.std(mat1[i,:])*len(mat1[i,:]))
				normat2 = (mat2[i,:]-np.mean(mat2[i,:]))/(np.std(mat2[i,:]))
				crosscor[i,:] = np.correlate(normat1,normat2,'same')
		elif axis == 1:
			for i in range(l):
				normat1 = (mat1[:,i]-np.mean(mat1[:,i]))/(np.std(mat1[:,i])*len(mat1[:,i]))
				normat2 = (mat2[:,i]-np.mean(mat2[:,i]))/(np.std(mat2[:,i]))
				crosscor[:,i] = np.correlate(normat1,normat2,'same')
		dumpfiles(crosscor,name)
		return crosscor

## get the cross correlation between two matrices
def getCrossCorrelation(sig1,sig2,name='power_crosscor'):
	if sig1.shape[0] != sig2.shape[0]:
		print('# ERROR # :: Cant calculate cross-correlation of unequal arrays')
		sys.exit()
	else:
		# setup cross correlation array
		crosscor = np.zeros(sig1.shape[0])
		crosscor = np.correlate(sig1,sig2,'same')
		dumpfiles(crosscor,name)
		return crosscor
		
## plot the cross correlation as a heatmap, from getCrossCorrelation
def plotCrossCorrelationMat(crosscor,ylabel='y',xlabel='x'):
	fig,ax=plt.subplots(figsize=(7,5))
	im = ax.imshow(crosscor,**kwargs,cmap='bwr') # blue white red map to show -1, 0, 1
	ax.set_ylabel(ylabel,**tnrfont) ; plt.xlabel(xlabel,**tnrfont)
	plt.colorbar(im)
#	fig.savefig('CrossCorrelation.png')
	return fig,ax

def plotCrossCorrelation(xarr,crosscor,ylabel='y',xlabel='x'):
	fig,ax=plt.subplots(figsize=(6,4))
	ax.plot(xarr,crosscor)
	ax.set_ylabel(ylabel,**tnrfont) ; plt.xlabel(xlabel,**tnrfont)
	return fig,ax

### example with ntot and Ep
## get sim
#simloc = getSimulation('/storage/space2/phrmsf/lowres_D_He3/0_34_p_90')
#ind_lst = list_sdf(simloc)
#
## quantities and ion species
#name = 'density_energy'
#species = getIonSpecies(sdfread(0))
#quantities = ['Derived_Number_Density_'+i for i in species]
#masses = [getMass(i) for i in species]
#
#try: # load cross cor
#	crosscor = read_pkl(name)	
#except: # calculate cross cor
#	# get number density and energy
#	narr = []
#	try:
#		for j in range(len(species)):
#			narr.append(load_batch_fieldmatrix([],quantities[j],para=False))
#	except:
#		get_batch_fieldmatrix(ind_lst,quantities=quantities,load=False)
#		for j in range(len(species)):
#			narr.append(load_batch_fieldmatrix([],quantities[j],para=False))
#	
#	# weighted number density according to mass
#	ntot = np.mean([masses[i]*narr[i] for i in range(len(species))],axis=0) # mean across each position
#	
#Ep = read_pkl('Protons_KEdensmatrix')
#dEp = Ep-np.mean(Ep[0:10,:])
#
## calculate cross correlation
#crosscorr = getCrossCorrelation(ntot,dEp)
#plotCrossCorrelation(crosscor,ylabel='t',xlabel='x')