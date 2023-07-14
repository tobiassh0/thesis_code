from list_new import *
from matplotlib.pyplot import cm
import my_constants as const

kwargs = {'interpolation':'nearest','origin':'lower','aspect':'auto'}

C3s = ['2','5','10','15','40','64']
quantity = 'Magnetic_Field_Bz'
maj_species = 'Deuterons'
in_klimprime = 100
in_wlimprime = 50

#### first run analysis
#for C3 in C3s:
#	print('C3 :: '+C3)
#	sim_loc = getSimulation('/storage/space2/phrmsf/C3_comp/C3_'+C3)
#	ind_lst = list_sdf(sim_loc)
#	d = sdfread(0)
#	species = getIonSpecies(d)
#	va = getAlfvenVel(d)
#	wcmaj = getCyclotronFreq(d,maj_species)
#	Nx=len(getGrid(d)[0])
#	xi1, xi2, xi3 = getConcentrationRatios(d)
#	C0 = int(C3)/xi3
#	C1 = xi1*C0 ; C2 = xi2*C0
#	Carr = np.array([int(C1),int(C2),int(C3)])
#	## field
#	try:
#		times = read_pkl('times')
#		fm = load_batch_fieldmatrix(ind_lst,quantity)
#	except:
#		times, fm = get_batch_fieldmatrix(ind_lst,quantity=quantity,load=True)
#		dumpfiles(fm,'fieldmatrix_'+quantity)
#	fm = fm - np.mean(fm[0:10,:])	# delta Bz

#	## energy
#	energy_int=energies(sim_loc,frac=1,plot=True,integ=True)

#	# fv_vA
#	_ = fv_vA(sim_loc,species_lst=species,para=False)

#	## FFT
#	try:
#		FT2d = read_pkl('FT_2d_'+quantity)
#	except:
#		FT2d = get2dTransform(fm,window=True)
#		dumpfiles(FT2d,'FT_2d_'+quantity)

#	klim = 0.5*2*const.PI*len(getGrid(sdfread(0))[0])/getGridlen(sdfread(0))
#	wlim = 0.5*2*const.PI*len(times)/times[-1]
#	klim_prime = klim*va/wcmaj
#	wlim_prime = wlim/wcmaj
#	w_lim, k_lim = FT2d.shape[0]*(in_wlimprime/wlim_prime), FT2d.shape[1]*(in_klimprime/klim_prime)

#	_,_ = power(klim_prime=klim_prime,wlim_prime=wlim_prime,wmax=20,kmax=40,wnorm=wcmaj,\
#			norm_omega=getOmegaLabel('Deuterons'),quantity='Magnetic_Field_Bz',plot=True)

#	FT2d = FT2d[:int(w_lim),:int(k_lim)]
#	fig,ax=plot2dTransform(FT2d, va_wci=[va,True], klim=in_klimprime, wlim=in_wlimprime, Omega_label=getOmegaLabel(maj_species),cmap='magma')
#	fig.savefig('FT2d_'+quantity+'.png',bbox_inches='tight')
##	im = axs[c].imshow(np.log10(FT2d),**kwargs,extent=[0,in_klimprime,0,in_wlimprime])


############################################################################################################


### compare all power spectra
c=0
fig,ax=plt.subplots(figsize=(10,4))
colors = cm.viridis(np.linspace(0.1, 1, len(C3s))) # can use any heatmap name here... 'cm.jet', 'cm.Accent', 'cm.bone_r' etc.
for C3 in C3s:
	print('C3 :: '+C3)
	sim_loc = getSimulation('/storage/space2/phrmsf/C3_comp/C3_'+C3)
	ind_lst = list_sdf(sim_loc)
	d = sdfread(0)
	species = getIonSpecies(d)
	va = getAlfvenVel(d)
	wcmaj = getCyclotronFreq(d,maj_species)
	Nx=len(getGrid(d)[0])
	xi1, xi2, xi3 = getConcentrationRatios(d)
	C0 = int(C3)/xi3
	C1 = xi1*C0 ; C2 = xi2*C0
	Carr = np.array([int(C1),int(C2),int(C3)])
	## times
	try:
		times = read_pkl('times')
	except:
		times, _ = get_batch_fieldmatrix(ind_lst,quantity=quantity,load=True)

	klim = 0.5*2*const.PI*len(getGrid(sdfread(0))[0])/getGridlen(sdfread(0))
	wlim = 0.5*2*const.PI*len(times)/times[-1]
	klim_prime = klim*va/wcmaj
	wlim_prime = wlim/wcmaj

	try:
		omegas_power = read_pkl('omegas_power') ; power = 10**read_pkl('log10_power')
	except:
		omegas_power, log10_power = power(klim_prime=klim_prime,wlim_prime=wlim_prime,wmax=20,kmax=40,wnorm=wcmaj,\
				norm_omega=getOmegaLabel('Deuterons'),quantity='Magnetic_Field_Bz',plot=True)
		power = 10**log10_power

	ax.plot(omegas_power/wcmaj,power,color=colors[c],label=C3)
	c+=1

for i in range(21):
	ax.axvline(i,linestyle='--',color='darkgrey')
ax.set_xlabel(r'$\omega/\Omega_D$',fontsize=20)
ax.set_ylabel('Power',fontsize=20)
plt.legend(loc='upper center',ncol=len(C3s))
ax.set_yscale('log')
fig.savefig('/storage/space2/phrmsf/C3_comp/Power_C3_comp.png',bbox_inches='tight')
#plt.show()
############################################################################################################


#### compare all fv_vA
#c=0
#width = 17 ; height = width/const.g_ratio
#fig,axs=plt.subplots(figsize=(width,height),nrows=3,ncols=len(C3s),sharex=True,sharey='row')
#fig.subplots_adjust(hspace=0.1,wspace=0.1)
#im = axs.copy()
#ylims=[(0.056,0.067),(0.044,0.054),(0.,0.6)]
#clims=[(2.8,4.0),(2.8,4.0),(1.0,4.5)]
#for C3 in C3s:
#	print('C3 :: '+C3)
#	sim_loc = getSimulation('/storage/space2/phrmsf/C3_comp/C3_'+C3)
#	ind_lst = list_sdf(sim_loc)
#	d = sdfread(0)
#	species = getIonSpecies(d)
#	va = getAlfvenVel(d)
#	wcmaj = getCyclotronFreq(d,maj_species)
#	Nx=len(getGrid(d)[0])
#	xi1, xi2, xi3 = getConcentrationRatios(d)
#	C0 = int(C3)/xi3
#	C1 = xi1*C0 ; C2 = xi2*C0
#	Carr = np.array([int(C1),int(C2),int(C3)])
#	## field
#	try:
#		times = read_pkl('times')
#	except:
#		times, _ = get_batch_fieldmatrix(ind_lst,quantity=quantity,load=True)

##	# fv_vA
##	_ = fv_vA(sim_loc,species_lst=species,para=False)
#	i = 0
#	for spec in species:
#		fv_spec = read_pkl('fv_'+spec)
#		v_spec  = read_pkl('v_'+spec)
#		vmin = np.min(v_spec) ; vmax = np.max(v_spec)
#		im[i,c]=axs[i,c].imshow(np.log10(fv_spec.T),**kwargs,extent=[0,times[-1]*wcmaj/(2*const.PI),vmin/va,vmax/va],cmap='jet',clim=clims[i])
#		axs[i,c].text(0.95,0.9,spec+' : '+str(Carr[i]),ha='right',va='top',transform=axs[i,c].transAxes)
#		axs[i,c].set_ylim(ylims[i])

#	c+=1

#axs[1,0].set_ylabel(r'$v/v_A$',fontsize=20)
#fig.text(0.5,0.035,r'$t\Omega_D/2\pi$',fontsize=20)
#p0 = axs[0,len(C3s)-1].get_position().get_points().flatten()
#p1 = axs[1,len(C3s)-1].get_position().get_points().flatten()
#p2 = axs[2,len(C3s)-1].get_position().get_points().flatten()
#print(p0,p1,p2)
#ax0_cbar = fig.add_axes([0.91, p1[1], 0.005 , p0[3]-p1[1]]) # [left bottom width height]
#ax1_cbar = fig.add_axes([0.91, p2[1], 0.005 , p2[3]-0.1]) # [left bottom width height]
#plt.colorbar(im[0,2], cax=ax0_cbar, orientation='vertical')
#plt.colorbar(im[2,2], cax=ax1_cbar, orientation='vertical')
#fig.text(0.95,(p0[3]+p1[1]-0.1)/2,r'$\log[f_i(v)]$',rotation='vertical',fontsize=20)
#fig.text(0.95,(p2[3])/2,r'$\log[f_\alpha(v)]$',rotation='vertical',fontsize=20)
##plt.show()
#fig.savefig('/storage/space2/phrmsf/C3_comp/fv_vA_collection.png',bbox_inches='tight')

