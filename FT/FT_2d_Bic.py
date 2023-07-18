
from func_load import *

sim_lst = ['traceT_0_00','traceT_D_99_T_01','traceT_D_89_T_11','traceT_0_50']
MuArr = [0,0.01,0.11,0.5]
Lambda = 2e-3
fig, axs = plt.subplots(figsize=(15,10),nrows=2,ncols=len(sim_lst),sharex='row',sharey='row')
fig.subplots_adjust(hspace=.45,wspace=0.1)
i=0
for sim in (sim_lst):
	## setup
	simLoc = getSimulation('/storage/space2/phrmsf/'+sim)
	print(sim)
	FT_2d = read_pkl('FT_2d_Magnetic_Field_Bz')
	times = read_pkl('times')
	bic = read_pkl('Bicohmat_ka_wa_0.06_40_nfft_2400_noverlap_1200')
	dt = times[-1]/len(times)
	dx = getdxyz(sdfread(0))
	LDe = getDebyeLength(sdfread(0),'Electrons')
	## freq limits
	klim = 2*0.5*const.PI/dx
	wlim = 2*0.5*const.PI/dt
	vA = getAlfvenVel(sdfread(0))
	wcyc = getCyclotronFreq(sdfread(0),'Deuterons')
	wnorm = wcyc
	knorm = 1/LDe#wcyc/vA
	klim_prime = klim/knorm
	wlim_prime = wlim/wnorm
	(nw,nk) = FT_2d.shape
	kmax = 0.06; wmax = 50
	print('kmax : {} , wmax :{}'.format(kmax,wmax))
	## freqs
	maj_species, maj2_species, min_species = getIonSpecies(sdfread(0))
	wce  = getCyclotronFreq(sdfread(0),'Electrons',Z=1)
	wpe	 = getPlasmaFreq(sdfread(0),species='Electrons')
	wpi  = getPlasmaFreq(sdfread(0),species='Deuterons')
	## LH
	Mu = MuArr[i]
	WLH = LowerHybridMassEffective(sdfread(0),(1-Mu-2*Lambda),Mu,'Deuterons','Tritons',wpe,wce,getMeanquantity(sdfread(0),'Derived_Number_Density_Electrons')) # effective
	wlh = np.sqrt(((wpi**2)+(wcyc**2))/(1+((wpe**2)/(wce**2))))
	axs[0,i].axhline(WLH/wnorm,linestyle='--',color='k',alpha=1.)
	## W2
	w2 = np.sqrt(wce*wcyc*((wce*wcyc + wpe**2)/((wce**2)+(wpe**2))))
	axs[0,i].axhline(w2/wnorm,linestyle='-.',color='k',alpha=1.)
	## chopping and plotting
	FT_2d = FT_2d[:int(nw*wmax/wlim_prime),:int(nk*kmax/klim_prime)]
	_,_,bic,_ = plot_bicoh(bic,extent=[0,kmax,0,kmax],bispectrum=False,smooth=True,cbar=False,clim=(0,1),cmap='jet')
	im0=axs[0,i].imshow(np.log10(FT_2d),interpolation='nearest',cmap='magma',origin='lower',aspect='auto',extent=[0,kmax,0,wmax],vmin=-4,vmax=6)
	for k in range(bic.shape[0]):
		for l in range(bic.shape[1]):
			if not bool(bic[k,l]):
				bic[k,l] = -1
	print(bic)
	bic = np.ma.masked_where(bic < 0, bic)
	cmap = plt.cm.jet
	cmap.set_bad(color='white') # removes white area in Bicoherence
	im1=axs[1,i].imshow(bic,interpolation='nearest',cmap='jet',origin='lower',aspect='auto',extent=[0,kmax,0,kmax],vmin=0,vmax=1)
	axs[0,i].set_ylim(0,wmax) ;	axs[1,i].set_ylim(0,kmax)
	axs[0,i].set_xlim(0,kmax) ; axs[1,i].set_xlim(0,kmax)
	if i == 0:
		axs[0,i].set_ylabel(r'$\omega/\Omega_D$',fontsize=20)
		axs[1,i].set_ylabel(r'$k_2\lambda_{De}$',fontsize=20)
#	if i == 2:
#		fig.colorbar(im0,ax=axs[0,i],orientation="horizontal")
#		fig.colorbar(im1,ax=axs[1,i],orientation="horizontal")
	labels = [0.00,0.01,0.02,0.03,0.04,0.05,0.06]
	axs[0,i].set_xticklabels(labels,rotation=45)
	axs[1,i].set_xticklabels(labels,rotation=45)
	i+=1
	
p00 = axs[0,0].get_position().get_points().flatten()
p01 = axs[0,1].get_position().get_points().flatten() 
p02 = axs[0,2].get_position().get_points().flatten()
p03 = axs[0,3].get_position().get_points().flatten()
p10 = axs[1,0].get_position().get_points().flatten()
p11 = axs[1,1].get_position().get_points().flatten() 
p12 = axs[1,2].get_position().get_points().flatten()
p13 = axs[1,3].get_position().get_points().flatten()

#ax0_cbar = fig.add_axes([p00[0], 0.94, p02[2]-p00[0], 0.01]) # [left bottom width height]
#plt.colorbar(im0, cax=ax0_cbar, orientation='horizontal')
#ax1_cbar = fig.add_axes([p10[0], 0.47, p12[2]-p10[0], 0.01]) # [left bottom width height]
#plt.colorbar(im1, cax=ax1_cbar, orientation='horizontal')
ax0_cbar = fig.add_axes([p00[0], 0.94, p03[2]-p00[0], 0.01]) # [left bottom width height]
plt.colorbar(im0, cax=ax0_cbar, orientation='horizontal')
ax1_cbar = fig.add_axes([p10[0], 0.47, p13[2]-p10[0], 0.01]) # [left bottom width height]
plt.colorbar(im1, cax=ax1_cbar, orientation='horizontal')


fig.text(0.5, 0.51, r'$k\lambda_{De}$', va='center', rotation='horizontal',fontsize=20) # xlabel, top row
fig.text(0.5, 0.04, r'$k_1\lambda_{De}$', va='center', rotation='horizontal',fontsize=20) # xlabel, bottom row

#axs[1,1].set_xlabel(r'$k_1\lambda_{De}$',fontsize=20) # v_A/\Omega_D
#axs[0,1].set_xlabel(r'$k\lambda_{De}$',fontsize=20)
os.chdir('/storage/space2/phrmsf/paper/')
plt.show()
fig.savefig('FT2d_Bic.png')
fig.savefig('FT2d_Bic.eps')










