
def getSpatialCorrelation(d0,F1='Magnetic_Field_Bz',F2='Electric_Field_Ex',\
								norm=[None,None],phase_res=500,plot=True):
	"""
		In:
			d0 : zeroth file
			F1 & F2 : names of fields you want to correlate, default Bz and Ex
			norm : normalisations of those fields (square rooted) default Bz and Ex
			phase_res : resolution (No. points) to take spatial cross-correlation
			plot : boolean of whether to plot data or not
		Out:
			dphase : array of phase shifts (rads) over which Rtdx is calculated
			Rtdx : Spatial cross correlation matrix
	"""
	# default normalisations for default fields
	if norm == [None,None]:
		norm[0] = (const.mult_magn_field)**.5
		norm[1] = (const.mult_elec_field)**.5
		
	# load fm and normalise
	fm1 = load_batch_fieldmatrix([],F1)*norm[0]
	fm2 = load_batch_fieldmatrix([],F2)*norm[1]
	
	L = getGridlen(d0)
	dx = getdxyz(d0)
	times = read_pkl('times')
	species = getAllSpecies(d0)#getIonSpecies(sdfread(0))
	wcmaj = getCyclotronFreq(d0,species[0])
	vA = getAlfvenVel(d0)
	kmax = 0.5*2*const.PI/dx

	FT1d = read_pkl('FT_1d_Magnetic_Field_Bz')
	print(FT1d.shape)
	# extract out maximum wavenumber for all times
	karg = np.argmax(FT1d[1:,:],axis=1) 
	karr = np.linspace(0,kmax,FT1d.shape[1])
	del FT1d
	k_star = karr[karg]
	thresh = None #k_star < 200*wcmaj/vA
	k_star = stats.mode(k_star)[0][0] # k_star[thresh]
	dphase = np.linspace(-2*const.PI,2*const.PI,phase_res)
	# shift between +-2pi, in terms of physical spacing
	Dx = dphase/k_star
	
	try:
		Rtdx = read_pkl('Rtdx_{}_{}'.format(F1,F2))
	except:
		Rtdx = np.zeros((fm1.shape[0],len(Dx))) # time-space
		print('There are {} rolls to do.'.format(len(Dx)))	
		for i in range(len(Dx)):
			print('roll :: ',i)
			fmroll = np.roll(fm1,int(Dx[i]/dx),axis=1) # dont need to roll through whole domain 
			Rtdx[:,i] = np.sum(fmroll*fm2*dx,axis=1)
		dumpfiles(Rtdx,'Rtdx_{}_{}'.format(F1,F2))
	if plot:
		plotSpatialCorrelation(F1,F2,Rtdx,dphase,times,species[-1])

	return dphase, Rtdx


def plotSpatialCorrelation(F1,F2,Rtdx,dphase,times,minspec):
	print(Rtdx.shape)
	T = times[-1]	
	tcmin = 2*const.PI/getCyclotronFreq(sdfread(0),minspec)
	fig,ax=plt.subplots(figsize=(8,6))
	ax.axvline(-const.PI,color='darkgrey',linestyle='--')
	ax.axvline(const.PI,color='darkgrey',linestyle='--')
	ax.imshow((Rtdx),**kwargs,cmap='jet',extent=[dphase[0],dphase[-1],0,T/tcmin])
	ax.set_xlabel('Phase shift, 'r'$\Delta \phi$'+'  '+'[rad]',**tnrfont)
	ax.set_xticks([-2*const.PI,-const.PI,0,const.PI,2*const.PI])
	ax.set_xticklabels([r'$-2\pi$',r'$-\pi$',r'$0$',r'$\pi$',r'$2\pi$'])
	ax.set_ylabel(r'$t$'+getOmegaLabel(minspec)+r'$/2\pi$',**tnrfont)
	fig.savefig('Rtdx_{}_{}.png'.format(F1[-2:],F2[-2:]),bbox_inches='tight')
	return None

#simloc=getSimulation('/storage/space2/phrmsf/ECRH/ECRH_JT60U_5_2')

if __name__=='__main__':
	from func_load import *
	os.chdir('/storage/space2/phrmsf/lowres_D_He3/')
	home = os.getcwd()
	sims = np.sort([i for i in os.listdir() if 'p_90' in i])
	hlabels = [int(i[2:4]) for i in sims]
	for sim in sims:
		simloc = getSimulation(sim)
		times = read_pkl('times')
		T = times[-1]
		tcmin = 2*const.PI/getCyclotronFreq(sdfread(0),'Protons')
		fields = ['Magnetic_Field_Bz','Electric_Field_Ex']#getFields()
		for i in range(len(fields)):
			for j in range(i,len(fields)): # one less field than previous loop (excludes repeats)
				F1 = fields[i]
				F2 = fields[j]
				_,_,norm1=Endict(F1)
				_,_,norm2=Endict(F2)
				dphase,Rtdx=getSpatialCorrelation(d0=sdfread(0),F1=F1,F2=F2,norm=[norm1,norm2],plot=True)
		os.chdir(home)