from func_load import *

def paraVelocity(INDEX):
	index, minority, mMin = INDEX
	print(index)
	return np.sqrt(2*getQuantity1d(sdfread(index),'Derived_Average_Particle_Energy_'+minority)/mMin)

def fv_vA(sim_loc,minority='Alphas',nval=10000,para=True):
	ind_lst = list_sdf(sim_loc)
	time = read_pkl('times')
	
	Nx = len(getGrid(sdfread(0))[0])
	Nt = len(time)
	vA = getAlfvenVel(sdfread(0))
	mMin = getMass(minority)
	vMin = np.zeros((Nt,Nx))

	## load velocity of Alphas
	try:
		vMin = read_pkl('v_'+minority)
	except:
		tind_lst = np.zeros((len(ind_lst),3),dtype='object')
		for i in range(len(ind_lst)):
			tind_lst[i,:] = [ind_lst[i], minority, mMin]
		if para:
			pool = mp.Pool(mp.cpu_count()//2)
			vMin = np.vstack(np.array(pool.map_async(paraVelocity,tind_lst).get(99999)))
			pool.close()
		else:
			vMin = np.zeros((Nt,Nx))
			for t in range(Nt):
				if 100*t/Nt%5==0:print(100*t/Nt, ' %')
				vMin[t,:] = np.sqrt(2*getQuantity1d(sdfread(t),'Derived_Average_Particle_Energy_'+minority)/mMin)
		dumpfiles(vMin,'v_'+minority)

	L = getGridlen(sdfread(0))
	tcMin = 2*const.PI/getCyclotronFreq(sdfread(0),minority)
	T = time[-1]/tcMin
	vmin  = np.min(vMin/vA) ; vmax = np.max(vMin/vA) 
	try:
		fMin = read_pkl('fv_'+minority)
	except:
		vMin = vMin/vA
		print('vmin ',vmin,'vmax ',vmax)
#		vMin = np.linspace(vmin,vmax,nval)
		fMin = np.zeros((Nt,nval))
		for t in range(fMin.shape[0]):
			xarr = np.linspace(0,L,Nx)
			yarr = vMin[t,:]
			fMin[t,:],_,_ = np.histogram2d(xarr,yarr,range=[[0,L],[vmin,vmax]],bins=(1,nval),density=True)
		#	yarr = np.linspace(vmin,vmax,nval)
		#	plt.plot(yarr,fAlpha[t,:]) ; plt.show()

		dumpfiles(fMin,'fv_'+minority)

	fig,ax = plt.subplots(figsize=(8,8/const.g_ratio))
	im = ax.imshow(np.log10(fMin.T),**kwargs,extent=[0,T,vmin,vmax],cmap='jet')

	cbar = plt.colorbar(im)
	cbar.ax.set_ylabel(r'$\log_{10}(f_\alpha)$', rotation=90)
	ax.set_xlabel(r'$t$'+getOmegaLabel(minority)+r'$/2\pi$',fontsize=14)	
#	ax.set_xlabel(r'$t/\tau_{cmin}$',fontsize=14)
	ax.set_ylabel(r'$v/v_A$',fontsize=14)
	plt.gca().ticklabel_format(useOffset=False)
	fig.savefig('f'+minority+'_v_vA.png',bbox_inches='tight')
#	plt.show()
	#plt.imshow((vMin/vA),**kwargs,cmap='seismic',extent=[0,L,0,T])
	#plt.colorbar()
	#plt.show()
	return fMin

#C3s = ['2','5']
#for C3 in C3s:
#	print('C3_'+C3)
#	sim_loc = getSimulation('/storage/space2/phrmsf/C3_comp/C3_'+C3)
#	fAlpha = fv_vA(sim_loc)

sim_loc = getSimulation('/storage/space2/phrmsf/old/JET_26148')
fAlpha = fv_vA(sim_loc,'Alphas')



