
# compares the power spectral density of multiple sims, recalcualtes power spectra to make consistent # TODO change this
def power_compare(sims,labels=[None],normspecies='Deuterons',wkmax=[25,40],quantity='Magnetic_Field_Bz',colors=None):
	fig,ax = plt.subplots(figsize=(10,5))
	if not colors: # check if self-defined colors
		colors = plt.cm.rainbow(np.linspace(0,1,len(sims)))
	#alphas = np.linspace(0.25,1.0,len(sims))	
	
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
		
		# re-calculate power spectra
		FT_2d = read_pkl('FT_2d_'+quantity)
		log10_power,omegas=powerspectrum(FT_2d,wnorm,[wlim_prime,klim_prime],[0,wkmax[0],0,wkmax[1]])
		power = 10**log10_power
		dw = (omegas[-1]-omegas[0])/len(omegas)
		
		# threshold		
		thresh = omegas < wnorm*wkmax[0]
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
	ax.set_xlim(0,wkmax[0])
	ax.set_xticks(np.linspace(0,wkmax[0],11))
	ax.legend(loc='best')
	fig.savefig('power_compare.png',bbox_inches='tight')
	plt.show()
	return None


if __name__=='__main__':
	from func_load import *
	os.chdir('/storage/space2/phrmsf/lowres_D_He3/')
	sims = np.sort([i for i in os.listdir() if 'p_90' in i])
	hlabels = np.array([int(i[2:4]) for i in sims])
	power_compare(sims,labels=hlabels,wkmax=[20,45],normspecies='Protons')
