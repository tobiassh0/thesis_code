
# Fourier transform of dBz in k-space at a given (normalised) time
def ENERGY_K(sim):
	sim_loc = getSimulation(sim)
	#ind_lst = sdf_list(sim_loc)
	
	quant = 'Magnetic_Field_Bz'
	#FT = read_pkl('FT_1d_'+quant)
	fm = load_batch_fieldmatrix([],'fieldmatrix_'+quant)
	times = read_pkl('times')
	(nt,nx) = fm.shape
	print(nt,nx)
	
	mean_Bz = np.mean(fm[0:100,:])
	print('fm mean Bz',mean_Bz)
	fm = ((fm-mean_Bz)**2)/(2*const.mu0)
	
	species = getIonSpecies(sdfread(0))
	tc = 2*const.PI/(getCyclotronFreq(sdfread(0),species[0]))
	L = getGridlen(sdfread(0))
	dx = L/nx
	wcyc = getCyclotronFreq(sdfread(0),species[0])
	va = getAlfvenVel(sdfread(0))
	knorm = wcyc/va
	klim = 0.5*2*const.PI*(1/dx)
	print('wcyc/va :: ',knorm)
	
	## Calculate 1d FT
	print('calculating 1d FT...')
	try:
		FT = read_pkl('FT_1d_edens_k')
	except:
		preshiftFT = np.zeros((nt, nx),complex)
		for t in range(0, nt):
			preshiftFT[t][:]=np.fft.fft(fm[t][:])
		FT = np.abs(np.fft.fftshift(preshiftFT,1))
		(nt,nk) = FT.shape
		FT = FT[:,int(FT.shape[1]/2):] ## all time but k > 0
		dumpfiles(FT,'FT_1d_edens_k')
		(nt,nk) = FT.shape
		del preshiftFT
	del fm 
	(nt,nk) = FT.shape
	print('nt :: {}, nk :: {}'.format(nt,nk))
	
	## extract time (tind)
	tn = 3.5 # normalised
	tind = nt*(tn/(times[-1]/tc))
	k = np.linspace(0,klim/knorm,nx//2)
	FT=(FT[int(tind),:]).flatten()
	
	## removing under-lying 1/k noise
	print('fitting noise...')
	#def noise(k,a):
	#	return a/k
	#popt, pcov = curve_fit(noise, k, FT)#, bounds=(1E6,1E7))
	#print(*popt)
	Enoise  = 10**5.5/k
	FT = FT-Enoise
	fig, ax1 = plt.subplots(figsize=(10,8))
	plt.plot(k,FT,color='b',label='data')
	ax1.set_ylabel(r'$(\Delta B_z)^2/2\mu_0$',fontsize=18)
	ax1.set_xlabel(r'$kv_A/\Omega_D$',fontsize=18)	
	ax1.set_xlim(0,120)
	ax1.set_ylim(-20000,150000)
	#ax1.set_yscale('log')
	ax1.set_title(r'$10^{5.5}/k$'+'  noise removed',fontsize=18)
	#plt.legend(loc='best')
	left, bottom, width, height = [30/140,100000/150000,50/140,30000/150000]
	ax2 = fig.add_axes([left,bottom,width,height])
	thresh = k < 40
	ax2.plot(k[thresh],FT[thresh])
	ax2.set_xlim(0,40)
	ax2.set_ylim(-20000,40000)
	plt.show()

if __name__=='__main__':
	from func_load import *
	ENERGY_K(sim='/storage/space2/phrmsf/old/resolve_e_gyro')
