
from func_load import *

def getSpatialCorrelation(d0,F1='Magnetic_Field_Bz',F2='Electric_Field_Ex',norm=[np.sqrt(const.mult_magn_field),np.sqrt(const.mult_elec_field)],plot=True):
	# in:
		# d0 : zeroth file
		# F1 & F2 : names of fields you want to correlate
		# norm : normalisations of those fields 	
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
	karg = np.argmax(FT1d[1:,:],axis=1) # extracts out maximum wavenumber for all times
	karr = np.linspace(0,kmax,FT1d.shape[1])
	del FT1d
	k_star = karr[karg]
	thresh = None #k_star < 200*wcmaj/vA
	k_star = stats.mode(k_star)[0][0] # k_star[thresh]
	dphase = np.linspace(-2*const.PI,2*const.PI,500)
	Dx = dphase/k_star # shift between +-2pi, in terms of physical spacing
	
	try:
		Rtdx = read_pkl('Rtdx_{}_{}'.format(F1,F2))
	except:
		Rtdx = np.zeros((fm1.shape[0],len(Dx))) # time-space
		for i in range(len(Dx)):
			print('roll :: ',i)
			fmroll = np.roll(fm1,int(Dx[i]/dx),axis=1)
			Rtdx[:,i] = np.sum(fmroll*fm2*dx,axis=1)
		dumpfiles(Rtdx,'Rtdx_{}_{}'.format(F1,F2))
	
	if plot:
		plotSpatialCorrealtion(F1,F2,Rtdx,dphase,times,species[-1])
	return dphase, Rtdx

def plotSpatialCorrealtion(F1,F2,Rtdx,dphase,times,minspec):
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
	ax.set_ylabel(r'$t$'+getOmegaLabel(minspec)+r'$2\pi$',**tnrfont)
	fig.savefig('Rtdx_{}_{}.png'.format(F1[-2:],F2[-2:]),bbox_inches='tight')
	return None

simloc=getSimulation('/storage/space2/phrmsf/ECRH/ECRH_JT60U_5_2')
fields = getFields()
for i in range(len(fields)):
	for j in range(i,len(fields)): # one less field than previous loop (excludes repeats)
		F1 = fields[i]
		F2 = fields[j]
		_,_,norm1=Endict(F1)	
		_,_,norm2=Endict(F2)
		dphase,Rtdx=getSpatialCorrelation(d0=sdfread(0),F1=F1,F2=F2,norm=[np.sqrt(norm1),np.sqrt(norm2)],plot=True)

