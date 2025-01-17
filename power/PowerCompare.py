
# compares the power spectral density of multiple sims, recalcualtes power spectra to make consistent # TODO change this
def power_compare(sims,labels=[None],normspecies='Deuterons',wkmax=[25,40],quantity='Magnetic_Field_Bz',colors=None,\
						xlims=[0,None],ylims=[None,None],leg=True,width=10,height=5,freqlabel=False,omegalabel=True):
	"""
		Compare the power spectra densities (PSD) of multiple sims through frequencies normalised to the normspecies.
		
		In:
			sims		 : the list of simulations for power spectra you want to compare
			labels		 : labels of each sim to put in the legend (see "leg" flag)
			normspecies	 : species which you want to normalise frequencies to
			wkmax	  	 : frequency and wavenumber maxima you want to calc the power spectra over [wmax, kmax]
			quantity	 : FT_2d quantity to load (defaults to dBz --> see batch_load for FT_2d calc)
			colors		 : color array of sims, defaults to plt.cm.rainbow 
			x&ylims		 : limits in x and y to plot, allows for a "zoomed in" version of the same spectra
			leg			 : boolean, legend on/off
			width/height : width and height of the figure (useful for zoomed plots)
			freqlabel	 : (default False) boolean to determine whether to plot as f/f_cminspec
			omegalabel	 : (default True) boolean to determine whether to plot as \omega/\Omega_minspec 
		Out: 
			Plots and saves the PSD comparison (on log scale). Returns "None"	
	"""	
	if freqlabel: 
		omegalabel = False # if boolean freqlabel on then turn off omegalabel
	# error check; opposite boolean flags
	if freqlabel == omegalabel:
		print('# ERROR # : trying to set ylabel as omega and freq')
		raise SystemExit

	wmax,kmax = wkmax
	xmin,xmax = xlims
	ymin,ymax = ylims
		
	## plot setup
	fig,ax = plt.subplots(figsize=(width,height))
	# colors
	if not colors: # check if self-defined colors
		colors = plt.cm.rainbow(np.linspace(0,1,len(sims)))
	# transparency
	# alphas = np.linspace(0.25,1.0,len(sims))	
	
	# maj/min harmonics
	if xmax != None:
		for j in range(xmin,xmax,1):
			ax.axvline(j,color='darkgrey',linestyle='--') # color='orange'
	else:
		for j in range(0,wmax,1):
			ax.axvline(j,color='darkgrey',linestyle='--')

	home=os.getcwd()
	for i in range(len(sims)):
	# for i in [0,2,4,6]:
		# load sim location		
		sim_loc = getSimulation(sims[i])
		times = read_pkl('times')	
		d0 = sdfread(1) # 0
		
		# frequency space limits and normalisation		
		wnorm = getCyclotronFreq(d0,normspecies)
		knorm = wnorm/getAlfvenVel(d0)
		klim_prime = (0.5*2*const.PI/getdxyz(d0))/knorm		
		wlim_prime = (0.5*2*const.PI/getdt(times))/wnorm
		
		# re-calculate power spectra for ranges given
		FT_2d = read_pkl('FT_2d_'+quantity)
		log10_power,omegas=powerspectrum(FT_2d,wnorm,[wlim_prime,klim_prime],[0,wmax,0,kmax])
		power = 10**log10_power
		dw = (omegas[-1]-omegas[0])/len(omegas)
				
		# plot ## psd = power/dw
		ax.plot(omegas/wnorm,np.log10(power/dw),label=labels[i],color=colors[i])#,alpha=alphas[i])
		os.chdir(home)	

	del times
	
	# formatting
	ax.set_ylabel('PSD',**tnrfont) # log-scale so unitless
	
	# frequency labels
	if omegalabel:
		ax.set_xlabel(r'$\omega/$'+getOmegaLabel(normspecies),**tnrfont)	
	if freqlabel:
		ax.set_xlabel(r'$f/$'+getFreqLabel(normspecies),**tnrfont)	
	# ax.set_yscale('log')
	
	# legend
	if leg:
		ax.legend(loc='upper left',labelspacing=0.1,borderpad=0.1,columnspacing=0.1,ncol=len(sims))
	
	# axis limits
	if ymin != None: # y-limits
		ax.set_ylim(ymin,ymax)
	if xmax == None: # x-limits
		print('# xmax is None #')
		ax.set_xlim(0,wmax)
		lims = (0,wmax)
	else:
		ax.set_xlim(xmin,xmax)
		lims = (xmin,xmax)
	
	# ax.set_xticks(np.linspace(0,wmax,6))
	ax.locator_params(axis='x',nbins=6)
	plt.show()
	fig.savefig('power_compare_{}_{}_short.png'.format(lims[0],lims[1]),bbox_inches='tight')	

	return None


if __name__=='__main__':
	from func_load import *

	## D-T
	# os.chdir('/storage/space2/phrmsf/traceT/')
	# # sims = np.sort([i for i in os.listdir() if 'traceT' in i])
	# # # rearrange so biggest to smallest T-concentration
	# # sims = np.append(sims,sims[0])
	# # sims = sims[1:]
	# # hlabels = np.array([int(str(i[-2] + i[-1])) for i in sims])
	# #
	# sims = ['traceT_D_50_T_50','traceT_D_89_T_11','traceT_D_99_T_01','traceT_D_100_T_00','cold_JET26148']
	# hlabels = [r'$50\%$',r'$11\%$',r'$1\%$',r'$0\%$','Baseline']
	# power_compare(sims,labels=hlabels,wkmax=[25,45],normspecies='Alphas',colors=['b','g','r','darkturquoise','k'],\
	# 				xlims=[10,25],leg=False,height=3)
	# sys.exit()

	## D-He3
	os.chdir('/storage/space2/phrmsf/lowres_D_He3/')
	# os.chdir('/run/media/phrmsf/My Passport/simulations/D-He3/pklfiles-lowres_D_He3/')
	sims = np.sort([i for i in os.listdir() if 'p_90' in i])
	sims = sims[1:] # remove 0%
	sims = sims[::2]
	hlabels = np.array([int(i[2:4]) for i in sims])	
	print(sims, hlabels)
	power_compare(sims,labels=hlabels,wkmax=[21,45],normspecies='Protons',xlims=[12,20],\
					omegalabel=True,leg=True,height=3)
