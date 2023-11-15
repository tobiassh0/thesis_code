
def func_exp(x,a,x0,sigma,offset):
	return a*np.exp(-(x-x0)**2/(2*sigma**2))+offset

def pearsonscc(sims,xiHe3):
	home = os.getcwd()
	# compare power spectra trend
		# shared area offset, relevant to 5% concentration
			# control
	wcp = 3.7*const.qe/getMass('Protons') # Bq/m
	simloc = getSimulation(sims[0])
	power5 = 10**read_pkl('log10_power')
	omegas5 = read_pkl('omegas_power')
	dw = omegas5[-1]/len(omegas5) 	
	power5/=dw # PSD
	thresh = omegas5/wcp < 20
	power5 = power5[thresh]
	omegas5= omegas5[thresh]
	#plt.plot(omegas5/wcp,np.log10(power5),label=hlabel[0])
	os.chdir(home)
	peaks=[0]
	colors=['k']
	perrs=[0]
	
	for i in range(1,len(sims)):
		print(sims[i])
		simloc = getSimulation(sims[i])
		power = 10**read_pkl('log10_power')
		omegas = read_pkl('omegas_power')
		dw = omegas[-1]/len(omegas) 	
		power/=dw # PSD
		thresh = omegas/wcp < 20
		power = power[thresh]
		omegas= omegas[thresh]
		w = np.linspace(-omegas[-1]/2,omegas[-1]/2,len(omegas))/wcp
		# get cross correlation
		crosscor = cc.getCrossCorrelation(np.log10(power),np.log10(power5)) # allows for correct difference between curves
		# weighted arithmetic mean
		mean = np.sum(w*crosscor)/np.sum(crosscor) ; sigma = 5#np.sqrt(np.sum(crosscor*(w-mean)**2)/np.sum(crosscor))
		offset = -10
		# fit an exp curve (normalised to wcp)
		popt, pcov = curve_fit(func_exp,w,crosscor,p0=[max(crosscor),mean,sigma,offset],maxfev=10000) # exponential fitting
		print(*popt)
		# append fit parameters and errors
		peaks.append(popt[1])# a, x0, sigma, offset
		perr = np.sqrt(np.diag(pcov))# one std on parameters
		perrs.append(perr[1])
		colors.append('b')
	#	plt.plot(w,crosscor)
	#	plt.plot(w,func_exp(w,*popt))
		os.chdir(home)
	
	fig,ax=plt.subplots(figsize=(7,5))
	for c,_x,_y,_yerr in zip(colors,xiHe3,peaks,perrs):
		if _x == xiHe3[0]:
			_yerr = None
		ax.errorbar(_x,_y,yerr=_yerr,xerr=None,fmt='o',color=c,capsize=15)
	r = np.corrcoef(xiHe3,peaks)[0,1]
	ax.set_ylabel(r'$\delta\omega/\Omega_p$',**tnrfont)
	ax.set_xlabel(r'$\xi_{He3}$',**tnrfont)
	ax.annotate(r'$r={}$'.format(np.around(r,2)),xy=(0.03,0.9),xycoords='axes fraction',ha='left',va='bottom',fontsize=18)
	ax.set_xlim(0.,0.5) ; ax.set_ylim(-0.05,0.45)
	xiHe3_err = np.std(xiHe3)
	perrs[0] = np.std(peaks) # re-write error on 5% xiHe3
	# fit ODR line of best fit with errors
	params, params_err = ODR_fit(x=xiHe3,y =peaks,sx=xiHe3_err,sy=perrs,curve='linear') # check list_new
	x_fit = np.linspace(0,.5,100)
	y_fit = func_linear(params,x_fit)
	ax.plot(x_fit,y_fit,color='r',linestyle='--')
	ax.annotate(r'$d(\delta \omega)/d \xi_{He3}=$'+r'$({}\pm{})$'.format(np.around(params[0],2),np.around(params_err[0],2))+r'$\Omega_p$',\
					xy=(0.95,0.05),xycoords='axes fraction',ha='right',va='bottom',fontsize=18)
	fig.savefig('PearsonsCorrCoef.png',bbox_inches='tight')
	plt.show()
	#np.savetxt('pearsonscrosscor.txt',np.array([xiHe3,peaks,perrs]).T)

def nEcc(sims,xiHe3):
	# cross correlation between number density and energy of protons
	for i in range(len(sims)):
		print(sims[i])
		simloc = getSimulation(sims[i])
		ind_lst = list_sdf(simloc)
		
		# quantities and ion species
		name = 'density_energy'
		species = getIonSpecies(sdfread(0))
		quantities = ['Derived_Number_Density_'+i for i in species]
		masses = [getMass(i) for i in species]
		
		try: # load cross cor
			raise SystemError
			crosscorr = read_pkl(name)	
		except: # calculate cross cor
			# get number density and energy
			narr = []
			try:
				for j in range(len(species)):
					narr.append(load_batch_fieldmatrix([],quantities[j],para=False))
			except:
				get_batch_fieldmatrix(ind_lst,quantities=quantities,load=False)
				for j in range(len(species)):
					narr.append(load_batch_fieldmatrix([],quantities[j],para=False))
			
			# weighted number density according to mass
			ntot = np.mean([masses[i]*narr[i] for i in range(len(species))],axis=0) # mean across each position and time
			plt.imshow(ntot-np.mean(ntot[0:10,:]),**kwargs,cmap='jet')
			plt.xlabel('x',**tnrfont) ; plt.xlabel('t',**tnrfont)
			plt.colorbar()
			plt.savefig(home+'/NumberDensity_{}.png'.format(xiHe3[i]))
			plt.clf()
			Ep = read_pkl('Protons_KEdensmatrix')
			dEp = Ep-np.mean(Ep[0:10,:])
			dEp = dEp/np.max(np.abs(dEp))
			# calculate cross correlation
			crosscorr = getCrossCorrelationMat(ntot,dEp)
		# plot cross correlation
		fig,ax=plotCrossCorrelationMat(crosscorr,ylabel='t',xlabel='x')
		os.chdir(home)
		fig.savefig('CrossCorrelation_{}.png'.format(xiHe3[i]))
		plt.clf()


if __name__ == '__main__':
	from func_load import *
#	from scipy.odr import *
	import correlation.cross_correlation as cc
	os.chdir('/storage/space2/phrmsf/lowres_D_He3/')
	sims = np.sort([i for i in os.listdir() if 'p_90' in i])
	xiHe3 = np.array([int(i[2:4])/100 for i in sims])
	pearsonscc(sims,xiHe3)
	#nEcc(sims,xiHe3)
	
	