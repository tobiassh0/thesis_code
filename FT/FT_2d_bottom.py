from list_new import *
import my_constants as const

sim_lst = ['traceT_0_00','traceT_D_99_T_01','traceT_D_89_T_11','traceT_0_50']
MuArr = [0,0.01,0.11,0.5]
Lambda = 2e-3

fig, axs = plt.subplots(figsize=(15,5),nrows=1,ncols=len(sim_lst),sharex='row',sharey='row')
fig.subplots_adjust(hspace=.45,wspace=0.1)
i=0

for sim in (sim_lst):
	## setup
	simLoc = getSimulation('/storage/space2/phrmsf/'+sim)
	print(sim)
	times = read_pkl('times')
#	bic = read_pkl('Bicohmat_ka_wa_0.06_40_nfft_2400_noverlap_1200')
	bic = read_pkl('Bicohmat_ka_wa_60_30_nfft_2400_noverlap_1200')
	dt = times[-1]/len(times)
	dx = getdxyz(sdfread(0))
	LDe = getDebyeLength(sdfread(0),'Electrons')
	## freq limits
	klim = 2*0.5*const.PI/dx
	wlim = 2*0.5*const.PI/dt
	vA = getAlfvenVel(sdfread(0))
	wcyc = getCyclotronFreq(sdfread(0),'Deuterons')
	wnorm = wcyc
	knorm = wcyc/vA#1/LDe
	klim_prime = klim/knorm
	wlim_prime = wlim/wnorm
	kmax = 60; wmax = 30 #0.025; 40
	print('kmax : {} , wmax :{}'.format(kmax,wmax))
	## freqs
	maj_species, maj2_species, min_species = getIonSpecies(sdfread(0))
	wce  = getCyclotronFreq(sdfread(0),'Electrons',Z=1)
	wpe	 = getPlasmaFreq(sdfread(0),species='Electrons')
	wpi  = getPlasmaFreq(sdfread(0),species='Deuterons')
	
	_,_,bic,_ = plot_bicoh(bic,extent=[0,kmax,0,kmax],bispectrum=False,smooth=True,cbar=False,clim=(0,1),cmap='jet')
	for k in range(bic.shape[0]):
		for l in range(bic.shape[1]):
			if not bool(bic[k,l]):
				bic[k,l] = -1
	bic = np.ma.masked_where(bic < 0, bic)
	cmap = plt.cm.jet
	cmap.set_bad(color='white') # removes white area in Bicoherence
	im1=axs[i].imshow(bic,interpolation='nearest',cmap='jet',origin='lower',aspect='auto',extent=[0,kmax,0,kmax],vmin=0,vmax=1) #[0,0.06,0,0.06]
	axs[i].set_ylim(0,40)
	axs[i].set_xlim(0,40)
	if i == 0:
		axs[i].set_ylabel(r'$k_2v_A/\Omega_D$',fontsize=20) # \lambda_{De}
#	labels = [0,0.005,0.01,0.015,0.02,0.025]
	labels = [0,5,10,15,20,25,30,35,40]
	axs[i].set_xticks(labels)
	axs[i].set_xticklabels(labels,rotation=45)
	axs[i].set_yticks(labels)
	axs[i].set_yticklabels(labels)
	axs[i].tick_params(axis='both',direction='out',top=False,right=False,left=True,bottom=True)
	i+=1
	
p0 = axs[0].get_position().get_points().flatten()
p1 = axs[1].get_position().get_points().flatten() 
p2 = axs[2].get_position().get_points().flatten()
p3 = axs[3].get_position().get_points().flatten()

ax1_cbar = fig.add_axes([p0[0], 0.97, p3[2]-p0[0], 0.02]) # [left bottom width height]
plt.colorbar(im1, cax=ax1_cbar, orientation='horizontal')

#fig.text(0.5, -0.03, r'$k_1\lambda_{De}$', va='center', rotation='horizontal',fontsize=20) # xlabel, bottom row
fig.text(0.5, -0.03, r'$k_1v_A/\Omega_D$', va='center', rotation='horizontal',fontsize=20) # xlabel, bottom row

#axs[1,1].set_xlabel(r'$k_1\lambda_{De}$',fontsize=20) # v_A/\Omega_D
#axs[0,1].set_xlabel(r'$k\lambda_{De}$',fontsize=20)
os.chdir('/storage/space2/phrmsf/paper/remake/')
#plt.show()
fig.savefig('FT2d_bottom.png',bbox_inches='tight')
fig.savefig('FT2d_bottom.eps',bbox_inches='tight')


