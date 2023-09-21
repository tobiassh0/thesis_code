## python files
import my_constants as const
from list_new import *

class Simulation():
	def __init__(self):
		import ASCII_logo
		del ASCII_logo
#		self.sim_file_loc = getSimulation('') # allows user to input the file destination in the dir where batch is run
		self.sim_file_loc = getSimulation('/storage/space2/phrmsf/lowres_D_He3/0_34_p_90')
		self.quantity = 'Magnetic_Field_Bz'
#		self.quantity = 'Derived_Number_Density_Deuterons'
		self.index_list = list_sdf(self.sim_file_loc)
		print(str(len(self.index_list))+' files')

#		times = np.zeros(len(self.index_list))
#		for i in range(len(self.index_list)):
#			file = sdfread(i)
#			times[i] = file.__dict__['Header']['time']
#		dumpfiles('times',times)

		try: ## load 0th and last file for later use in calculation of constants
			print('Loading first and last file...')
			self.file0 = sdfread(0)
			self.filelast = sdfread(self.index_list[-1])
			print('Done.')
		except:
			## loads 00001 if 00000 is not present
			print('# ERROR # : Couldn\'t load files.')
			raise SystemExit

		try:
			self.times = read_pkl('times')
		except:
			print('# ERROR # : Couldn\'t load times.')

		# grid_length, N_grid_points, duration, N_time_points = batch_getSimRanges((self.index_list,self.file0,self.filelast,self.times))
		self.Nx = len(getGrid(self.file0)[0]) # can generalise to use length of array (i.e. Ex, Bz or number dens)

		## Creates list of available fieldmatrix pickled files in the sim directory
		lst_dir = os.listdir(self.sim_file_loc)
		pkl_list = [item for item in lst_dir if self.quantity+'.pkl' in item] # check if fieldmatrix of quantity is pickled
		lst_FT2d = [item for item in pkl_list if 'FT_2d' in item]
		lst_Fm = [item for item in pkl_list if 'fieldmatrix' in item]
		self.quantities = getFields()
		
		## calc and load fm & FT2d
		if len(pkl_list) == 0: # dont have pkl files # should be on first running
			self.times, self.fieldmatrix = get_batch_fieldmatrix(self.index_list,self.quantities,load=True)
			if 'times.pkl' not in os.listdir() or (self.times==None).any():
				dumpfiles(self.times, 'times')
		else: # have pkl files
			pkl_file = self.quantity
			Nos = ['no','No','n','','N']
			if len(lst_FT2d)>0:print('You have [{}] available FT 2d pkl file(s)'.format(len(lst_FT2d)))
			if len(lst_Fm)>0:print('You have [{}] available fieldmatrix pkl file(s)'.format(len(lst_Fm)))
			if pkl_file not in Nos:
				try:
					self.FT_2d = read_pkl('FT_2d_'+self.quantity)
					print('shape FT2d (w,k) :: {}'.format(np.shape(self.FT_2d)))
					self.fieldmatrix = load_batch_fieldmatrix(self.index_list,self.quantity,para=False) 
					load_fm = False # default dont load as already worked #  fm if loading FT2d
				except: 
					load_fm = True
				if load_fm: self.fieldmatrix = load_batch_fieldmatrix(self.index_list,self.quantity,para=False) # load fm as well
			else:
				raise SystemExit

		## delta Bz
		print(self.fieldmatrix.shape)
		if self.quantity == 'Magnetic_Field_Bz':
			self.fieldmatrix = self.fieldmatrix - np.mean(self.fieldmatrix[0:10,:]) # convert to delta Bz rather than pure Bz

		## Times and Ions
		self.T = self.times[-1]
		self.Nt = len(self.times)
		self.dt = (self.times[-1]-self.times[0])/self.Nt
		spec_names = getIonSpecies(self.file0) # get a list of species names from scanning the 0th file's keys for an always returned value "Number_Density_" + species
		print('# Ions :: ',len(spec_names), ' :: ', spec_names)
		maj_species, maj2_species, min_species = spec_names
		min_species = getAllSpecies(sdfread(0))[-1]
		Zmaj = getChargeNum(maj_species)
		Zmin = getChargeNum(min_species)

		## Find concentration ratios
		Xiarr = getConcentrationRatios(self.file0)
		Xi1, Xi2, Xi3 = Xiarr
		Xi2_Xi1 = Xi2/Xi1	;		Xi = Xi3/Xi1
		print('Xi1 :: {}\nXi2 :: {}\nXi2/Xi1 :: {}\nXi3 :: {}'.format(Xi1,Xi2,Xi2_Xi1,Xi3))
		conc_ratios = [Xi1,Xi2,Xi2_Xi1,Xi3,Xi]

		## Characteristic frequencies
		self.wc_maj  = getCyclotronFreq(self.file0,maj_species,Z=Zmaj)
		self.wc_min  = getCyclotronFreq(self.file0,min_species,Z=Zmin)
		self.wce  = getCyclotronFreq(self.file0,'Electrons',Z=1)
		self.va   = getAlfvenVel(self.file0)
		self.lambdaDe= getDebyeLength(self.file0,'Electrons')
		self.wpe	 = getPlasmaFreq(self.file0,species='Electrons')
		self.wpi  = getPlasmaFreq(self.file0,species=maj_species)
		self.tc_maj = 2*const.PI/self.wc_maj
		self.tc_min = 2*const.PI/self.wc_min

		## Lengths and Normalisation
		self.L = getGridlen(self.file0)
		self.rLe = getLarmorRadius(self.file0,'Electrons')
		self.rL_maj = getLarmorRadius(self.file0,maj_species)
		self.dx = getdxyz(self.file0)
		print('### RATIO rLe/dx :: {} \n rLmaj/dx :: {}'.format(self.rLe/self.dx,self.rL_maj/self.dx))
		self.tnorm = self.tc_min
		self.wnorm = self.wc_min
		self.knorm = self.wnorm/self.va #maj  #1/self.lambdaDe
		print('normalisation w: ', self.wnorm, ' [Hz]')
		
	### PLOT ENERGY DENSITIES OF FIELD COMPONENTS ### 
		energy_int=energies(sim_loc=self.sim_file_loc,frac=1,plot=True,integ=True)
#		energy_int = 0

	### PLOT CIGARETTE PLOTS ###
		ciggies(self.sim_file_loc,species_lst=getAllSpecies(self.file0),para=False) # doesnt return anything

	### FOURIER TRANSFORMS ###
		self.klim = 0.5*2*const.PI*self.Nx/self.L
		self.wlim = 0.5*2*const.PI*self.Nt/self.T
		self.klim_prime = self.klim/self.knorm
		self.wlim_prime = self.wlim/self.wnorm
		self.tlim_prime = self.T/self.tnorm

#		in_klimprime = float(input('!>> klim_prime maximum on plot (0 for max): ')) ## in normalised units
#		in_wlimprime = float(input('!>> wlim_prime maximum on plot (0 for max): ')) ## in normalised units
		in_klimprime = 40 #100 # use when lazy
		in_wlimprime = 25 #50
		if in_klimprime == 0:
			in_klimprime = self.klim_prime
		if in_wlimprime == 0:
			in_wlimprime = self.wlim_prime
		in_tlimprime = self.tlim_prime # maximum time

	### Saving Simulation Data as txt and pkl ### 
		simdatadict={'species' : spec_names,
				'Xi1, Xi2, Xi2_Xi1, Xi3, Xi' : conc_ratios,
				'va' : self.va,
				'wlim' : self.wlim,
				'klim' : self.klim,
				'wlim_prime' : self.wlim_prime,
				'klim_prime' : self.klim_prime,
				'wci' : self.wc_maj,
				'wpi' : self.wpi,
				'wce' : self.wce,
				'wpe' : self.wpe,
				'lambda_De' : self.lambdaDe,
				'tlim_prime' : self.tlim_prime ,
				'dt' : self.dt,
				'T' : self.T,
				'dx' : self.dx,
				'L' : self.L,
				'rLmaj/dx' : self.rL_maj/self.dx,
				'rLe/dx' : self.rLe/self.dx,
				'energy_integral' : energy_int,
		         }

		dumpfiles(simdatadict,'Sim_Data') # dump as pkl so can read directly later on
		###
		def dict_dump_txt(dict_data,filename='Sim_Data.txt'):
			print('Writing sim data...')
			with open(filename, 'w') as f: 
				for key, value in dict_data.items(): 
					f.write('%s:%s\n' % (key, value))
			print('Done.')
			return None
		###
		dict_dump_txt(simdatadict) # write as .txt so can physically read, since pkl writes in binary


	### FT 1d ###
		try:
			self.FT_1d = read_pkl('FT_1d_'+self.quantity)
			(nt,nk) = self.FT_1d.shape ; print('FT_1d shape (t,k):: {}'.format(nt,nk))
			self.FT_1d = self.FT_1d[:int(nt*in_tlimprime/self.tlim_prime),:int(nk*in_klimprime/self.klim_prime)]
		except:
			print('Creating 1d FFT...')
			self.FT_1d = get1dTransform(self.fieldmatrix,window=False)
			dumpfiles(self.FT_1d,'FT_1d_'+self.quantity)
			## Limits the size of the loaded array to our inputted limits so that it normalises properly ##
			k_lim, t_lim = self.FT_1d.shape[1]*(in_klimprime/self.klim_prime), self.FT_1d.shape[0]*(in_tlimprime/self.tlim_prime)
			print(k_lim,t_lim,self.FT_1d.shape)
			print('k_file_lim ',k_lim,'t_file_lim ', t_lim)
			self.FT_1d = self.FT_1d[:int(t_lim),:int(k_lim)]
		fig_1, ax_1 = plot1dTransform(self.FT_1d,klim=in_klimprime,tlim=in_tlimprime,wlabel=getOmegaLabel(min_species))
		plotting(fig_1,ax_1,'FT_1d_'+self.quantity)
		del self.FT_1d

	### FT 2d ###
		try:
			self.FT_2d = read_pkl('FT_2d_'+self.quantity)
			w_lim, k_lim = self.FT_2d.shape[0]*(in_wlimprime/self.wlim_prime), self.FT_2d.shape[1]*(in_klimprime/self.klim_prime)
			self.FT_2d = self.FT_2d[:int(w_lim),:int(k_lim)]
		except:
			print('Creating all FT_2ds...')
			for quant in self.quantities: # create all FT_2d arrays for available fields
				fmq = load_batch_fieldmatrix([],quant)
				FT_2d = get2dTransform(fmq,window=True)
				dumpfiles(FT_2d,'FT_2d_'+quant)
			# load the one you want to analyse
			self.FT_2d = read_pkl('FT_2d_'+self.quantity)
			# limits the size of the loaded array to our inputted limits so that it normalises properly ##
			w_lim, k_lim = self.FT_2d.shape[0]*(in_wlimprime/self.wlim_prime), self.FT_2d.shape[1]*(in_klimprime/self.klim_prime)
			self.FT_2d = self.FT_2d[:int(w_lim),:int(k_lim)]
		print('plotting shape: ',np.shape(self.FT_2d))
		fig, ax = plot2dTransform(self.FT_2d,klim=in_klimprime,wlim=in_wlimprime,klabel=getWavenumberLabel(min_species),wlabel=getOmegaLabel(min_species))#min_species

	### COLD PLASMA DISPERSION ###
#		Te = getTemperature('Electrons') ; Ti = getTemperature('Deuterons') 
#		k, omegas, w_LH = calcNewColdDisp(in_klimprime=100,Te=Te,Ti=Ti)
		omegas = self.wnorm*np.linspace(0,in_wlimprime,10000) # TODO: find out why doesn't go through 0 (Im and Re components)
		k1,k2,k3=coldplasmadispersion(self.file0,omegas=omegas) # two solutions to the cold plasma dispersion
		k1,k2,k3 = k1/self.knorm, k2/self.knorm, k3/self.knorm # can plot all three but not needed, k2 is main real branch
		ax = ColdWaveModes(ax,[self.wpi, self.wpe, self.wc_maj, self.wce],self.va,self.wnorm,LHact=True)#,W2=True) ## check list_new to see which modes we can plot (bool)
		## W_LH = LowerHybridMassEffective(self.file0,self.wpe,self.wce,n_e)
		thresh = k2 > 0 # threshold the array so it only plots the FAW and not the horizontal line to the 0 parts of the dispersion
		ax.plot(k2[thresh],omegas[thresh]/self.wnorm,color='k',linestyle='-',alpha=0.75)
		ax.set_ylim(0,in_wlimprime) # cold dispersion goes above this so this cuts it off
		ax.set_xlim(0,in_klimprime) # " "
		plotting(fig,ax,'FT_2d_'+self.quantity)

	### POWER SPECTRUM ###
		_,_ = power(klim_prime=self.klim_prime,wlim_prime=self.wlim_prime,wmax=35,kmax=in_klimprime,wnorm=self.wnorm,\
			norm_omega=getOmegaLabel(min_species),quantity=self.quantity,plot=True)
	
#	### Poynting ###
#		FTSmag = Poynting2dFT(self.times,self.Nt,self.Nx,in_klimprime=in_klimprime,plot=True)
		
	### BICOHERENCE SPECTRA ###
		karea = 30
		warea = 30
		nfft = self.Nt//5
		noverlap = nfft//2
		bispec=True #boolean, True as default will calculate both bicoh AND bispec
#		klabel = r'$\lambda_{De}$'
		klabel = r'$v_A/\Omega_p$'
		try:
			bicname = 'Bicohmat_ka_wa_{}_{}'.format(karea,warea)
			bisname = 'Bispecmat_ka_wa_{}_{}'.format(karea,warea)
			## bicoherence
			bic = read_pkl(bicname)
			fig,ax,_,_=plot_bicoh(bic,extent=[0,karea,0,karea],bispectrum=False,smooth=True,cbar=True,clim=(0,1),cmap='jet')
			ax.set_xlabel(r'$k_1$'+klabel,fontsize=18)
			ax.set_ylabel(r'$k_2$'+klabel,fontsize=18)
			plotting(fig,ax,'bicoherence_karea_{}_warea_{}'.format(karea,warea))
			## bispectrum
			bis = read_pkl(bisname)
			fig,ax,_,_=plot_bicoh(bis,extent=[0,karea,0,karea],bispectrum=True,smooth=True,cbar=True,cmap='jet')
			ax.set_xlabel(r'$k_1$'+klabel,fontsize=18)
			ax.set_ylabel(r'$k_2$'+klabel,fontsize=18)
			plotting(fig,ax,'bispectrum_karea_{}_warea_{}'.format(karea,warea))
		except:
			fig, ax = getBicoh(karea,warea,self.fieldmatrix,self.dt,self.T,self.L,self.wnorm,self.knorm,nfft=nfft,\
				noverlap=noverlap,window=True,bispectrum=bispec,klabel=klabel)


if __name__ == '__main__':
	Simulation()
