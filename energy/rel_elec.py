from func_load import *

simloc = getSimulation('/storage/space2/phrmsf/ECRH/ECRH_JT60U_5')
d0 = sdfread(0)
species = ['Electrons','FElectrons','Protons']
equant = [i+'_KEdens' for i in species]
times = read_pkl('times')
labels=[r'$e$',r'$e_{rel}$',r'$p$']
tcp = 2*const.PI/getCyclotronFreq(d0,'Protons')
wce = getCyclotronFreq(d0,'Electrons')
kce = wce/getAlfvenVel(d0)

### Energies
## load energies
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

## FT2d
FT2d = read_pkl('FT_2d_Magnetic_Field_Bz')
klim = 0.5*2*const.PI/getdxyz(d0)
wlim = 0.5*2*const.PI/getdt(times)
extents = [0,klim/kce,0,wlim/wce]
plt.imshow(np.log10(FT2d),**kwargs,extent=extents,clim=(-6,6))
plt.axhline(1,color='k',linestyle='--')
plt.show()
