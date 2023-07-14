from list_new import *
import my_constants as const
from scipy.special import gamma
import os, sys

kwargs = {'interpolation':'nearest','origin':'lower','aspect':'auto'}


sim_lst = ['traceT_0_50']#,'traceT_D_70_T_30','traceT_D_82_T_18','traceT_D_89_T_11','traceT_D_95_T_05','traceT_D_99_T_01','traceT_0_00']

for sim in sim_lst:
	sim_loc = getSimulation('/storage/space2/phrmsf/'+sim)
	species = getIonSpecies(sdfread(0))
	wcmaj = getCyclotronFreq(sdfread(0),species[0])
	omegas = wcmaj*np.linspace(0,20,10000)
	try:
		FT2d = read_pkl('FT2d_ReIm')
	except:
		fieldmatrix = load_batch_fieldmatrix()
		FT2d = get2dTransform(fieldmatrix)
		dumpfiles(FT2d,'FT2d_ReIm')	
	ReFT2d = np.real(FT2d)
	ImFT2d = np.imag(FT2d)
	try:
		GAMMA = read_pkl('GAMMA_Z')
	except:
		GAMMA = np.zeros((ReFT2d.shape[0],ReFT2d.shape[1]))
		for l in range(ReFT2d.shape[0]):
			for k in range(ReFT2d.shape[1]):
				GAMMA[l,k] = np.abs(gamma(ReFT2d[l,k]+1j*ImFT2d[l,k]))
		print(GAMMA.shape)
		dumpfiles(GAMMA,'GAMMA_Z')
	print(GAMMA)
	plt.imshow(np.log10(GAMMA),**kwargs)
	plt.show()


#	phase = np.arctan(ImFT2d/ReFT2d)
#	plt.imshow(phase,**kwargs,cmap='jet')	
#	plt.show()
	
#	fig = plt.figure(figsize=(10,10)) ; ax = fig.gca(projection='3d')
#	ax.plot_trisurf()
#	plt.show()
#	ax.plot_surface()
#	ax.set_xlabel('Re')
#	ax.set_ylabel('Im')
#	ax.set_zlabel('Phase')
#	plt.show()
