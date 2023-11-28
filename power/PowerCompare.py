
# compares the power spectral density of multiple sims, recalcualtes power spectra to make consistent # TODO change this
def power_compare(sims,labels=[None],normspecies='Deuterons',wkmax=[25,40],quantity='Magnetic_Field_Bz',colors=None,\
						xlims=[0,None],ylims=[None,None],leg=True,width=10,height=5):
	"""
		Compare the power spectra densities (PSD) of multiple sims through frequencies normalised to the normspecies.
		
		In:
			sims			 : the list of simulations for power spectra you want to compare
			labels		 : labels of each sim to put in the legend (see "leg" flag)
			normspecies	 : species which you want to normalise frequencies to
			wkmax			 : frequency and wavenumber maxima you want to calc the power spectra over
			quantity		 : FT_2d quantity to load (defaults to dBz --> see batch_load for FT_2d calc)
			colors		 : color array of sims, defaults to plt.cm.rainbow 
			x&ylims		 : limits in x and y to plot, allows for a "zoomed in" version of the same spectra
			leg			 : boolean, legend on/off
			width/height : width and height of the figure (useful for zoomed plots)
		Out: 
			Plots and saves the PSD comparison (on log scale). Returns "None"	
	"""	
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
	for i in range(0,wmax,1):
		ax.axvline(i,color='darkgrey',linestyle='--')

	i=0
	home=os.getcwd()	
	for sim in sims:
		# load sim location		
		sim_loc = getSimulation(sim)
		times = read_pkl('times')	
		d0 = sdfread(0)
		
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
		
		# threshold		
		thresh = omegas < wnorm*wmax
		omegas = omegas[thresh]/wnorm
		psd = power[thresh]/dw
		
		# plot
		ax.plot(omegas,psd,label=labels[i],color=colors[i])#,alpha=alphas[i])
		i+=1
		os.chdir(home)	

	del times
	# formatting
	ax.set_ylabel('PSD',**tnrfont) ; ax.set_xlabel(r'$\omega/$'+getOmegaLabel(normspecies),**tnrfont)	
	ax.set_yscale('log')
	ax.set_xticks(np.linspace(0,wmax,6))
	if leg: # legend
		ax.legend(loc='best',labelspacing=0.1,borderpad=0.1,ncol=1)
	# axis limits
	if ymin != None: # y-limits
		ax.set_ylim(ymin,ymax)
	if not xmin: # x-limits
		print(' xmin is None ')		
		ax.set_xlim(0,wmax)
		fig.savefig('power_compare_{}_{}.png'.format(xmin,wmax),bbox_inches='tight')	
	else:
		ax.set_xlim(xmin,xmax)
		fig.savefig('power_compare_{}_{}.png'.format(xmin,xmax),bbox_inches='tight')	
	
	#plt.show()
	return None


if __name__=='__main__':
	from func_load import *
	## D-T 
	os.chdir('/storage/space2/phrmsf/traceT/')
#	sims = np.sort([i for i in os.listdir() if 'traceT' in i])
#	# rearrange so biggest to smallest T-concentration
#	sims = np.append(sims,sims[0])
#	sims = sims[1:]
#	hlabels = np.array([int(str(i[-2] + i[-1])) for i in sims])
	sims = ['traceT_D_50_T_50','traceT_D_89_T_11','traceT_D_99_T_01','traceT_D_100_T_00','cold_JET26148']
	hlabels = [r'$50\%$',r'$11\%$',r'$1\%$',r'$0\%$','Cold']
	power_compare(sims,labels=hlabels,wkmax=[25,45],normspecies='Alphas',colors=['b','g','r','darkcyan','k'])#,xlims=[10,25],leg=False,height=3)
	sys.exit()

	## D-He3
	os.chdir('/storage/space2/phrmsf/lowres_D_He3/')
	sims = np.sort([i for i in os.listdir() if 'p_90' in i])
	hlabels = np.array([int(i[2:4]) for i in sims])	
	power_compare(sims,labels=hlabels,wkmax=[25,45],normspecies='Protons')
