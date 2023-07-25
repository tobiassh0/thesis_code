
from func_load import *


sim_lst = ['traceT_0_00','traceT_D_99_T_01','traceT_D_89_T_11','traceT_0_50']
Xi2 = [0,0.01,0.11,0.5]
Lambda = 2e-3

fig, axs = plt.subplots(figsize=(15,5),nrows=1,ncols=len(sim_lst),sharex='row',sharey='row')
fig.subplots_adjust(hspace=.45,wspace=0.1)
i=0
for sim in (sim_lst):
	## setup
	simLoc = getSimulation('/storage/space2/phrmsf/'+sim)
	print(sim)
	FT_2d = read_pkl('FT_2d_Magnetic_Field_Bz')
	times = read_pkl('times')
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
	(nw,nk) = FT_2d.shape
#	kmax = 0.025; wmax = 30
	kmax = 30; wmax = 30
	print('kmax : {} , wmax :{}'.format(kmax,wmax))
	## freqs
	maj_species, maj2_species, min_species = getIonSpecies(sdfread(0))
	wce  = getCyclotronFreq(sdfread(0),'Electrons',Z=1)
	wpe	 = getPlasmaFreq(sdfread(0),species='Electrons')
	wpi  = getPlasmaFreq(sdfread(0),species='Deuterons')
#	## LH
#	xi2 = Xi2[i]
#	WLH = LowerHybridMassEffective(sdfread(0),wpe,wce,getMeanquantity(sdfread(0),'Derived_Number_Density_Electrons'))
#	axs[i].axhline(WLH/wnorm,linestyle='--',color='k',alpha=1.)
#	## W2
#	w2 = np.sqrt(wce*wcyc*((wce*wcyc + wpe**2)/((wce**2)+(wpe**2))))
#	axs[i].axhline(w2/wnorm,linestyle='-.',color='k',alpha=1.)
	## chopping and plotting
	FT_2d = FT_2d[:int(nw*wmax/wlim_prime),:int(nk*kmax/klim_prime)]
	im0=axs[i].imshow(np.log10(FT_2d),interpolation='nearest',cmap='magma',origin='lower',aspect='auto',extent=[0,kmax,0,wmax],vmin=-2,vmax=6)
	axs[i].set_ylim(0,wmax)
	axs[i].set_xlim(0,kmax)
	if i == 0:
		axs[i].set_ylabel(r'$\omega/\Omega_D$',fontsize=20)
#	labels = [0.00,0.005,0.01,0.015,0.02,0.025]
	labels = [0,5,10,15,20,25,30]
	axs[i].set_xticklabels(labels,rotation=45)
	axs[i].tick_params(axis='both',direction='out',top=False,right=False,left=True,bottom=True)
	i+=1
	
p0 = axs[0].get_position().get_points().flatten()
p1 = axs[1].get_position().get_points().flatten() 
p2 = axs[2].get_position().get_points().flatten()
p3 = axs[3].get_position().get_points().flatten()
ax0_cbar = fig.add_axes([p0[0], 0.97, p3[2]-p0[0], 0.02]) # [left bottom width height]
plt.colorbar(im0, cax=ax0_cbar, orientation='horizontal')

#fig.text(0.5, -0.03, r'$k\lambda_{De}$', va='center', rotation='horizontal',fontsize=20) # xlabel, top row
fig.text(0.5, -0.03, r'$kv_A/\Omega_D$', va='center', rotation='horizontal',fontsize=20) # xlabel, top row
#fig.text(0.5, 0.04, r'$k_1\lambda_{De}$', va='center', rotation='horizontal',fontsize=20) # xlabel, bottom row

#axs[1,1].set_xlabel(r'$k_1\lambda_{De}$',fontsize=20) # v_A/\Omega_D
#axs[0,1].set_xlabel(r'$k\lambda_{De}$',fontsize=20)
os.chdir('/storage/space2/phrmsf/paper/remake/')
#plt.show()
fig.savefig('FT2d_top.png',bbox_inches='tight')
fig.savefig('FT2d_top.eps',bbox_inches='tight')
