
def CompareNxM_energy(sims):
	"""
	Plot a NxM panelled image of the energies of each sim in sims, with one label for each (if shared species)
		In:
			sims
		Out:
			
	"""
	# try 3, then 2 then 1x1 (i.e. single sim). If none of the above then opt for the least difference of 2 or 3
	try:
		N = np.array([3,2])
		M = np.array([3,2])
		n = len(sims) % N
		mind = np.argmin(n)
		print(M[mind], N[mind])
	except:
		N = M = 1 # single sim/image
	
	return None

if __name__=='__main__':
	from func_load import *
	os.chdir('/storage/space2/phrmsf/lowres_D_He3')
	sims = np.sort([i for i in os.listdir() if 'p_90' in i])
	labels = [i[2:4] for i in sims]
	CompareNxM_energy(sims)