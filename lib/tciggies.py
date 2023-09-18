from func_load import *

def ciggies_test(sim_loc,species_lst=['Deuterons','Alphas'],nval=10000,para=False,eload=True,vload=False):
	if eload:
		vload = False
	elif vload:
		eload = False
	elif eload == vload:
		print('## ERROR ## :: vload and eload are set equal') ; raise SystemExit
	else:
		print('## ERROR ## :: no value set for eload vload') ; raise SystemExit

	# initial parameters, dimensions and constants
	ind_lst = list_sdf(sim_loc)
	time = read_pkl('times')
	d0 = sdfread(0)
	Nx = len(getGrid(d0)[0])
	Nt = len(time)
	vA = getAlfvenVel(d0)
	species_lst = [i for i in species_lst if i != ''] # remove missing species
	
	# setup mass and velocity arrays
	for species in species_lst:
		L = getGridlen(d0)
		tcSpecies = 2*const.PI/getCyclotronFreq(d0,species)
		T = time[-1]/tcSpecies
		if vload: # velocities
			ylabel = r'$v/v_A$'
			cbar_label = r'$\log_{10}[f(v)]$'
			figname = 'fv_vA_'
			matnorm = vA
			try:
				fSpecies = read_pkl('fv_'+species)
				fload = True
			except:
				fload = False
		if eload: # energies
			ylabel = r'$E$' + '  ' + '['+r'$keV$'+']'
			cbar_label = r'$\log_{10}[f(E)]$'
			figname = 'fE_keV_'
			matnorm = 1000*const.qe
			try:
				fSpecies = read_pkl('fE_'+species)
				fload = True
			except:
				fload = False
		# check if fload is possible
		if fload:
			fSpecies = read_pkl(figname+species)
			if vload: # should load if here, can't be here otherwise
				matSpecies = read_pkl('v_'+species)
			if eload:
				matSpecies = read_pkl(species+'_KEdensmatrix')/dens
		else: # calculate distribution if can't load it already
			dens = getMeanquantity(sdfread(0),'Derived_Number_Density_'+species)
			massSpecies = getMass(species)
			matSpecies = np.zeros((Nt,Nx))
			try: # try loading
				if vload:
					matSpecies = read_pkl('v_'+species)
				if eload:
					matSpecies = read_pkl(species+'_KEdensmatrix')/dens
			except: # didn't load, calculate instead
				if eload:
					matSpecies,_ = getEnergies([species+'_KEdens'],[species],Nt,dump=True)/dens # loads energy density matrix of species
				if vload:
					tind_lst = np.zeros((len(ind_lst),3),dtype='object')
					for i in range(len(ind_lst)):
						tind_lst[i,:] = [ind_lst[i], species, massSpecies]
					if para: # parallel calculation, just for velocities ## TODO: change this to be more general
						pool = mp.Pool(mp.cpu_count()//2)
						matSpecies = np.vstack(np.array(pool.map_async(paraVelocity,tind_lst).get(99999)))
						pool.close()
					else:
						for t in range(Nt):
							if int(100*t/Nt)%5==0:print(100*t/Nt, ' %')
							matSpecies[t,:] = np.sqrt(2*getQuantity1d(sdfread(t),'Derived_Average_Particle_Energy_'+species)/massSpecies)		
					dumpfiles(matSpecies,'v_'+species)
			# initialise distribution matrix
			fSpecies = np.zeros((Nt,nval))
			# calculate fv/fE matrix
			for t in range(fSpecies.shape[0]): # loop through time
				xarr = np.linspace(0,L,Nx)
				yarr = matSpecies[t,:]
				fSpecies[t,:],_,_ = np.histogram2d(xarr,yarr,range=[[0,L],[matmin,matmax]],bins=(1,nval),density=True)
			# dump distribution matrix
			dumpfiles(fSpecies,figname+species)

		# normalise
		matSpecies = matSpecies/matnorm
		matmin  = np.min(matSpecies) ; matmax = np.max(matSpecies) ; matmean = np.mean(matSpecies)

		# extents of matrix per species
		extents = [0,T,matmin,matmax]
		print('min val ',matmin,'max mat ',matmax)

		# setup figure
		fig,ax = plt.subplots(figsize=(8,8/const.g_ratio))
		im = ax.imshow(np.log10(fSpecies.T),**kwargs,extent=extents,cmap='jet')

		# colorbar and labels
		cbar = plt.colorbar(im)
		cbar.ax.set_ylabel(cbar_label, rotation=90,fontsize=18)
		ax.set_xlabel(r'$t$'+getOmegaLabel(species)+r'$/2\pi$',fontsize=16)
		ax.set_ylabel(ylabel,fontsize=16)
		
		# set y-lim
		if (0.9*matmean) < matmin:
			if (1.1*matmean) > matmax:
				ax.set_ylim(matmin,matmax)
			else:
				ax.set_ylim(matmin,1.1*matmean)
		else:
			ax.set_ylim(0.9*matmean,1.1*matmean)

		plt.gca().ticklabel_format(useOffset=False)
		# savefig
		fig.savefig(figname+species+'.png',bbox_inches='tight')
	return None



sim_loc = getSimulation('/home/space/phrmsf/Documents/EPOCH/epoch-4.17.16/epoch1d/old/0006qp')
ciggies_test(sim_loc,['Deuterons','Electrons'])#,nval=10000,para=False,eload=True,vload=False)


