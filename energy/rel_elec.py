from func_load import *

simloc = getSimulation('/storage/space2/phrmsf/ECRH/ECRH_JT60U_5')
d0 = sdfread(0)
species = ['Electrons','FElectrons','Protons']
times = read_pkl('times')
labels=[r'$e$',r'$e_{rel}$',r'$p$']
wcp = getCyclotronFreq(d0,'Protons')
tcp = 2*const.PI/wcp
wce = getCyclotronFreq(d0,'Electrons')
tce = 2*const.PI/wce
tnorm = tcp
wnorm = 2*const.PI/tnorm
knorm = wnorm/getAlfvenVel(d0)
# freq limits
klim = 0.5*2*const.PI/getdxyz(d0)
wlim = 0.5*2*const.PI/getdt(times)

# ciggies 
ciggies(simloc,species)
sys.exit()
## energies
#equant = [i+'_KEdens' for i in species]
#_,_=getEnergies(equant,species,nt=len(times))
## plot
#for i in range(len(species)):
#	e_dens = read_pkl(equant[i])
#	plt.plot(times/tcp,e_dens-np.mean(e_dens[0:10]),label=labels[i])
#plt.legend(loc='best')
#plt.ylabel(r'$\Delta u $'+'  ['+r'$Jm^{-3}$'+']',**tnrfont)
#plt.xlabel(r'$t\Omega_p/2\pi$',**tnrfont)
#plt.xlim(0,times[-1]/tcp)
#plt.savefig('du_particles.png',bbox_inches='tight')
#plt.show()

## 2d FT
#fmEy = load_batch_fieldmatrix([],'Electric_Field_Ey')
#FT2d = get2dTransform(fmEy)
#dumpfiles(FT2d,'FT_2d_Electric_Field_Ey')
#FT2d = read_pkl('FT_2d_Electric_Field_Ey')
#extents = [0,klim/knorm,0,wlim/wnorm]
#plt.imshow(np.log10(FT2d),**kwargs,cmap='magma',extent=extents)#,clim=(-10,))
#plt.axhline(1,color='k',linestyle='--')
#plt.xlim(0,90)
#plt.ylim(0,100)
#plt.ylabel(r'$\omega/\Omega_p$',**tnrfont)
#plt.xlabel(r'$kv_A/\Omega_p$',**tnrfont)
#plt.savefig('FT2d_Bz_elec.png',bbox_inches='tight')
#plt.show()

# loop FT and power spectra
extents = [0,klim/knorm,0,wlim/wnorm]
quantities = getFields()
extents = [0,klim/knorm,0,times[-1]/tnorm]
for quant in quantities:
	fmEy = load_batch_fieldmatrix([],quant)
	FT1d = get1dTransform(fmEy)
	dumpfiles(FT1d,'FT_1d_'+quant)
	FT1d = read_pkl('FT_1d_'+quant)
	plt.imshow(np.log10(FT1d),**kwargs,cmap='bwr',extent=extents)#,clim=(-10,))
#	plt.axhline(1,color='k',linestyle='--')
	plt.xlim(0,90)
#	plt.xlim(0,0.05)
#	plt.ylim(0,100)
#	plt.ylabel(r'$\omega/\Omega_p$',**tnrfont)
	plt.ylabel(r'Time,  '+r'$\tau_{cp}$',**tnrfont)
	plt.xlabel(r'$kv_A/\Omega_p$',**tnrfont)
	plt.savefig('FT1d_'+quant+'.png',bbox_inches='tight')
#	power(klim_prime=klim/knorm,wlim_prime=wlim/wnorm,wmax=60,kmax=90,wnorm=wcp,norm_omega=r'$\Omega_p$',quantity=quant,plot=True,read=False)

